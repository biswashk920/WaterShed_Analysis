"""utils.py - shared helpers for the Kali watershed notebooks (02 onward)."""
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio import features
from rasterio.windows import Window
from shapely.geometry import shape, MultiPolygon, LineString
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, Rectangle
from matplotlib.ticker import FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

NODATA = -9999.0
plt.rcParams.update({"font.size": 8, "figure.dpi": 300})


# ---------------------------------------------------------------- paths
def get_paths():
    cwd = Path.cwd().resolve()
    root = cwd.parent if cwd.name == "Notebooks" else cwd
    p = dict(ROOT=root, INPUTS=root / "Inputs",
             RAS=root / "Outputs" / "Rasters", VEC=root / "Outputs" / "Vectors",
             FIG=root / "Outputs" / "Figures", TAB=root / "Outputs" / "Tables")
    for k in ("RAS", "VEC", "FIG", "TAB"):
        p[k].mkdir(parents=True, exist_ok=True)
    return SimpleNamespace(**p)


def wbt_run(label, func, *args, **kwargs):
    t0 = time.time()
    rc = func(*args, **kwargs)
    print(f"{label:<30} rc={rc}  {time.time() - t0:6.1f} s")
    assert rc == 0, f"{label} failed"


# ---------------------------------------------------------------- watershed grid
class WatershedGrid:
    """Bounding window of the watershed on the full cropped-DEM grid.
    read() returns any full-grid raster clipped to that window, NaN outside the basin."""

    def __init__(self, mask_path, boundary_path, outlet_path=None):
        with rasterio.open(mask_path) as src:
            m = src.read(1).astype(bool)
            self.full_transform = src.transform
            self.full_shape = src.shape
            self.crs = src.crs
            self.res = src.res
        rr = np.where(m.any(axis=1))[0]
        cc = np.where(m.any(axis=0))[0]
        self.window = Window(int(cc[0]), int(rr[0]),
                             int(cc[-1] - cc[0] + 1), int(rr[-1] - rr[0] + 1))
        self.transform = rasterio.windows.transform(self.window, self.full_transform)
        self.mask = m[rr[0]:rr[-1] + 1, cc[0]:cc[-1] + 1]
        del m
        h, w = self.mask.shape
        left, top = self.transform * (0, 0)
        right, bottom = self.transform * (w, h)
        self.extent = (left, right, bottom, top)
        self.cell_area_km2 = self.res[0] * self.res[1] / 1e6
        self.n_cells = int(self.mask.sum())
        self.area_km2 = self.n_cells * self.cell_area_km2
        self.boundary = gpd.read_file(boundary_path).to_crs(self.crs)
        self.outlet = None
        if outlet_path is not None:
            o = gpd.read_file(outlet_path).to_crs(self.crs)
            self.outlet = (o.geometry.iloc[0].x, o.geometry.iloc[0].y)

    def __repr__(self):
        return (f"WatershedGrid(window={self.mask.shape}, cells={self.n_cells:,}, "
                f"area={self.area_km2:,.1f} km2, res={self.res[0]:.2f} m)")

    def read(self, path, band=1):
        with rasterio.open(path) as src:
            same = (src.shape == self.full_shape and np.allclose(
                tuple(src.transform)[:6], tuple(self.full_transform)[:6], atol=1e-6))
            if not same:
                raise ValueError(f"{Path(path).name} is not on the dem_cropped grid")
            arr = src.read(band, window=self.window).astype("float32")
            nd = src.nodata
        if nd is not None:
            arr[arr == nd] = np.nan
        arr[~self.mask] = np.nan
        return arr

    def save(self, arr, path, nodata=NODATA, dtype="float32"):
        """Write a clipped array as a compressed GeoTIFF (NaN -> nodata)."""
        out = np.where(np.isfinite(arr), arr, nodata).astype(dtype)
        profile = dict(driver="GTiff", height=out.shape[0], width=out.shape[1], count=1,
                       dtype=dtype, crs=self.crs, transform=self.transform,
                       nodata=nodata, compress="lzw", tiled=True,
                       blockxsize=256, blockysize=256)
        with rasterio.open(path, "w", **profile) as dst:
            dst.write(out, 1)
        return path


# ---------------------------------------------------------------- classification
def classify(arr, breaks):
    """Inner break values -> class ids 1..len(breaks)+1 (NaN kept)."""
    out = np.full(arr.shape, np.nan, dtype="float32")
    ok = np.isfinite(arr)
    out[ok] = np.digitize(arr[ok], breaks) + 1
    return out


def class_table(cls, labels, ws):
    ids = np.arange(1, len(labels) + 1)
    cnt = np.array([(cls == i).sum() for i in ids])
    return pd.DataFrame({"class_id": ids, "label": labels, "cells": cnt,
                         "area_km2": (cnt * ws.cell_area_km2).round(2),
                         "percent": (100 * cnt / cnt.sum()).round(2)})


def sym_limits(arr, pct=98):
    v = float(np.nanpercentile(np.abs(arr), pct))
    return -v, v


def save_table(df, path, index=False):
    df.to_csv(path, index=index)
    return path


# ---------------------------------------------------------------- vectors
def classes_to_vector(cls, ws, labels, out_stem, sieve_px=50):
    """Classified raster -> polygons (one MultiPolygon per class). Saves .gpkg and .shp.
    Patches smaller than sieve_px cells are merged into neighbours first."""
    a = np.where(np.isfinite(cls), cls, 0).astype("int16")
    if sieve_px:
        a = features.sieve(a, size=sieve_px, mask=ws.mask, connectivity=4)
    recs = {}
    for geom, val in features.shapes(a, mask=ws.mask, transform=ws.transform, connectivity=4):
        recs.setdefault(int(val), []).append(shape(geom))
    ids = sorted(recs)
    gdf = gpd.GeoDataFrame({"class_id": ids,
                            "label": [labels[i - 1] for i in ids],
                            "geometry": [MultiPolygon(recs[i]) for i in ids]}, crs=ws.crs)
    gdf["area_km2"] = (gdf.area / 1e6).round(3)
    gdf.to_file(str(out_stem) + ".gpkg", driver="GPKG")
    gdf.to_file(str(out_stem) + ".shp")
    return gdf


def contours_to_vector(dem, ws, interval, out_stem, index_every=5):
    """Contour lines from a clipped DEM array. Saves .gpkg and .shp."""
    import contourpy
    x = ws.transform.c + (np.arange(dem.shape[1]) + 0.5) * ws.transform.a
    y = ws.transform.f + (np.arange(dem.shape[0]) + 0.5) * ws.transform.e
    gen = contourpy.contour_generator(x, y[::-1], np.ma.masked_invalid(dem[::-1]),
                                      line_type="Separate")
    lo = np.floor(np.nanmin(dem) / interval) * interval
    hi = np.ceil(np.nanmax(dem) / interval) * interval
    rows = []
    for lev in np.arange(lo, hi + interval, interval):
        for line in gen.lines(lev):
            if len(line) > 3:
                rows.append((float(lev), bool(round(lev / interval) % index_every == 0),
                             LineString(line)))
    gdf = gpd.GeoDataFrame({"elev": [r[0] for r in rows], "index": [r[1] for r in rows],
                            "geometry": [r[2] for r in rows]}, crs=ws.crs)
    gdf["geometry"] = gdf.geometry.simplify(ws.res[0] * 0.5)
    gdf.to_file(str(out_stem) + ".gpkg", driver="GPKG")
    gdf.to_file(str(out_stem) + ".shp")
    return gdf


# ---------------------------------------------------------------- map figure
def make_map(arr, ws, out_png, title, label, cmap="viridis", vmin=None, vmax=None,
             classes=None, hillshade=None, overlay=None, overlay_kw=None,
             show_outlet=True, dpi=300):
    """One standard map.
    Continuous : cmap, vmin/vmax (default 2-98 percentile); vertical colour bar outside, right.
    Categorical: classes=[(id, label, colour), ...] with consecutive ids 1..n; legend top-left.
    Always     : basin boundary, north arrow (top-right), scale bar (bottom-left), km axes."""
    step = max(1, int(np.ceil(max(arr.shape) / 3000)))
    a = arr[::step, ::step]
    left, right, bottom, top = ws.extent
    w_m, h_m = right - left, top - bottom
    fig_w = 8.0
    fig, ax = plt.subplots(figsize=(fig_w, float(np.clip(fig_w * h_m / w_m, 5, 11))))

    alpha = 1.0
    if hillshade is not None:
        ax.imshow(hillshade[::step, ::step], extent=ws.extent, cmap="gray",
                  interpolation="nearest")
        alpha = 0.6

    if classes:
        ids = [c[0] for c in classes]
        cmap_c = ListedColormap([c[2] for c in classes])
        norm = BoundaryNorm(np.array(ids + [ids[-1] + 1]) - 0.5, cmap_c.N)
        im = ax.imshow(a, extent=ws.extent, cmap=cmap_c, norm=norm,
                       interpolation="nearest", alpha=alpha)
    else:
        fin = a[np.isfinite(a)]
        lo = np.percentile(fin, 2) if vmin is None else vmin
        hi = np.percentile(fin, 98) if vmax is None else vmax
        im = ax.imshow(a, extent=ws.extent, cmap=cmap, vmin=lo, vmax=hi,
                       interpolation="nearest", alpha=alpha)

    ws.boundary.boundary.plot(ax=ax, color="k", linewidth=0.8, zorder=4)
    if overlay is not None:
        overlay.plot(ax=ax, **(overlay_kw or dict(color="k", linewidth=0.3)))
    if show_outlet and ws.outlet is not None:
        ax.plot(*ws.outlet, marker="*", ms=11, mfc="red", mec="k", mew=0.6, ls="none", zorder=6)
        ax.annotate("Outlet", ws.outlet, xytext=(7, 7), textcoords="offset points",
                    fontsize=7, fontweight="bold", zorder=6,
                    bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))

    # north arrow, top-right
    ax.annotate("", xy=(0.95, 0.945), xytext=(0.95, 0.865), xycoords="axes fraction",
                arrowprops=dict(facecolor="k", edgecolor="k", width=3, headwidth=9, headlength=9),
                zorder=9)
    ax.text(0.95, 0.952, "N", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=11, fontweight="bold", zorder=9)

    # categorical legend, top-left, small font
    if classes:
        handles = [Patch(fc=c[2], ec="0.3", lw=0.4, label=c[1]) for c in classes]
        leg = ax.legend(handles=handles, loc="upper left", fontsize=6, title=label,
                        title_fontsize=6.5, framealpha=0.92, borderpad=0.6,
                        handlelength=1.2, handleheight=0.9, labelspacing=0.35)
        leg.set_zorder(10)
    else:
        # vertical colour bar outside the map, right side, same height as the map
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="3%", pad=0.12)
        cb = fig.colorbar(im, cax=cax, orientation="vertical")
        cb.ax.tick_params(labelsize=6.5, length=2.5)
        cb.set_label(label, fontsize=7.5, labelpad=6)

    # scale bar, bottom-left
    nice = [1, 2, 5, 10, 20, 25, 50, 100, 200]
    km = max([n for n in nice if n <= w_m / 1000 * 0.2] or [1])
    f = km * 1000 / w_m
    x0, x1, yb = 0.04, 0.04 + f, 0.045
    ax.add_patch(Rectangle((x0 - 0.02, 0.015), f + 0.04, 0.07, transform=ax.transAxes,
                           fc="white", ec="0.4", lw=0.5, alpha=0.92, zorder=8))
    ax.plot([x0, x1], [yb, yb], transform=ax.transAxes, color="k", lw=2.2,
            solid_capstyle="butt", zorder=9)
    ax.text((x0 + x1) / 2, yb + 0.012, f"{km} km", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6.5, zorder=9)

    km_fmt = FuncFormatter(lambda v, _: f"{v / 1000:.0f}")
    ax.xaxis.set_major_formatter(km_fmt)
    ax.yaxis.set_major_formatter(km_fmt)
    ax.set_xlabel(f"Easting (km), EPSG:{ws.crs.to_epsg()}", fontsize=7)
    ax.set_ylabel("Northing (km)", fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)
    fig.savefig(out_png, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.show()
    plt.close(fig)
    return out_png

def map_furniture(ax, ws, title, outlet=True):
    """Boundary, outlet, north arrow (top-right), scale bar (bottom-left), km axes, title."""
    left, right, bottom, top = ws.extent
    w_m = right - left
    ws.boundary.boundary.plot(ax=ax, color="k", linewidth=0.8, zorder=4)
    if outlet and ws.outlet is not None:
        ax.plot(*ws.outlet, marker="*", ms=11, mfc="red", mec="k", mew=0.6, ls="none", zorder=6)
        ax.annotate("Outlet", ws.outlet, xytext=(7, 7), textcoords="offset points",
                    fontsize=7, fontweight="bold", zorder=6,
                    bbox=dict(fc="white", ec="none", alpha=0.7, pad=1))
    ax.annotate("", xy=(0.95, 0.945), xytext=(0.95, 0.865), xycoords="axes fraction",
                arrowprops=dict(facecolor="k", edgecolor="k", width=3, headwidth=9, headlength=9),
                zorder=9)
    ax.text(0.95, 0.952, "N", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=11, fontweight="bold", zorder=9)
    nice = [1, 2, 5, 10, 20, 25, 50, 100, 200]
    km = max([n for n in nice if n <= w_m / 1000 * 0.2] or [1])
    f = km * 1000 / w_m
    x0, x1, yb = 0.04, 0.04 + f, 0.045
    ax.add_patch(Rectangle((x0 - 0.02, 0.015), f + 0.04, 0.07, transform=ax.transAxes,
                           fc="white", ec="0.4", lw=0.5, alpha=0.92, zorder=8))
    ax.plot([x0, x1], [yb, yb], transform=ax.transAxes, color="k", lw=2.2,
            solid_capstyle="butt", zorder=9)
    ax.text((x0 + x1) / 2, yb + 0.012, f"{km} km", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=6.5, zorder=9)
    km_fmt = FuncFormatter(lambda v, _: f"{v / 1000:.0f}")
    ax.xaxis.set_major_formatter(km_fmt)
    ax.yaxis.set_major_formatter(km_fmt)
    ax.set_xlabel(f"Easting (km), EPSG:{ws.crs.to_epsg()}", fontsize=7)
    ax.set_ylabel("Northing (km)", fontsize=7)
    ax.tick_params(labelsize=6.5)
    ax.set_title(title, fontsize=11, fontweight="bold")
    ax.set_xlim(left, right)
    ax.set_ylim(bottom, top)