#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import LinearSegmentedColormap

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 20,
    "axes.labelsize": 20,
    "xtick.labelsize": 20,
    "ytick.labelsize": 20,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


def parse_args():
    p = argparse.ArgumentParser(description="Plot Hi-C triangle plus one-chromosome genomic feature tracks.")
    p.add_argument("--cool", required=True)
    p.add_argument("--chromosome", required=True)
    p.add_argument("--centromere", required=True)
    p.add_argument("--copia", required=True)
    p.add_argument("--gypsy", required=True)
    p.add_argument("--satdna", required=True)
    p.add_argument("--gene", required=True)
    p.add_argument("--expression", required=True)
    p.add_argument("--cg", required=True)
    p.add_argument("--chg", required=True)
    p.add_argument("--chh", required=True)
    p.add_argument("--out-prefix", required=True)
    p.add_argument("--gypsy-ltr")
    p.add_argument("--gypsy-int")
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--no-pdf", action="store_true", help="Do not save PDF.")
    p.add_argument("--hic-bin-step", type=int, default=1, help="Plot every Nth Hi-C bin for faster previews.")
    p.add_argument("--max-distance", type=float, default=None)
    p.add_argument("--hic-vmax-percentile", type=float, default=99.5)
    p.add_argument("--font-file", default=None, help="Optional Arial .ttf/.otf path. Use this to force true Arial.")
    return p.parse_args()


def red_cmap():
    return LinearSegmentedColormap.from_list(
        "white_to_red", [(0.0, (1, 1, 1)), (0.5, (1, 0, 0)), (1.0, (1, 0, 0))], N=256
    )


def read_track(path, chrom):
    df = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1, 2, 3], names=["chr", "start", "end", "value"])
    df = df[df["chr"] == chrom].copy()
    df[["start", "end", "value"]] = df[["start", "end", "value"]].apply(pd.to_numeric, errors="coerce")
    df = df.dropna()
    df["mid_mb"] = (df["start"] + df["end"]) / 2 / 1e6
    return df


def scale_mtg_tracks(copia, gypsy, satdna, gene, cg, chg, chh):
    # Keep the same scaling logic as MTG_02_plot_multi_tracks.R.
    copia = copia.copy()
    gypsy = gypsy.copy()
    satdna = satdna.copy()
    gene = gene.copy()
    cg = cg.copy()
    chg = chg.copy()
    chh = chh.copy()

    copia["value"] = np.minimum(copia["value"].to_numpy(dtype=float) * 1.8, 1.0)
    gypsy["value"] = np.minimum(gypsy["value"].to_numpy(dtype=float) * 1.5, 1.0)
    gene["value"] = np.minimum(gene["value"].to_numpy(dtype=float), 1.0)
    sat_max = float(np.nanmax(satdna["value"].to_numpy(dtype=float))) if not satdna.empty else 0.0
    if np.isfinite(sat_max) and sat_max > 0:
        satdna["value"] = satdna["value"].to_numpy(dtype=float) / sat_max
    cg["value"] = np.minimum(cg["value"].to_numpy(dtype=float), 1.0)
    chg["value"] = np.minimum(chg["value"].to_numpy(dtype=float), 1.0)
    chh["value"] = np.minimum(chh["value"].to_numpy(dtype=float), 1.0)
    return copia, gypsy, satdna, gene, cg, chg, chh


def read_expr(path, chrom, column):
    if path == "none":
        return pd.DataFrame(columns=["chr", "start", "end", "value", "mid_mb"])
    df = pd.read_csv(path, sep="\t")
    required = {"chr", "start", "end", column}
    if not required.issubset(df.columns):
        raise SystemExit(f"ERROR: expression file must contain columns: {', '.join(sorted(required))}")
    df = df[df["chr"] == chrom][["chr", "start", "end", column]].copy()
    df.columns = ["chr", "start", "end", "value"]
    df[["start", "end", "value"]] = df[["start", "end", "value"]].apply(pd.to_numeric, errors="coerce")
    df = df.dropna()
    df["mid_mb"] = (df["start"] + df["end"]) / 2 / 1e6
    return df


def read_ticks(path, chrom):
    if not path:
        return pd.DataFrame(columns=["chr", "start", "end", "mid_mb"])
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame(columns=["chr", "start", "end", "mid_mb"])
    df = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    df = df[df["chr"] == chrom].copy()
    df["mid_mb"] = (df["start"] + df["end"]) / 2 / 1e6
    return df


def read_centromere(path, chrom):
    df = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    df = df[df["chr"] == chrom]
    if df.empty:
        return None
    return float(df.iloc[0]["start"]) / 1e6, float(df.iloc[0]["end"]) / 1e6


def draw_track(ax, df, color, label, ymax, cen=None, fill=True):
    if cen:
        ax.axvline(cen[0], color="black", lw=1.0, ls=(0, (4, 4)))
        ax.axvline(cen[1], color="black", lw=1.0, ls=(0, (4, 4)))
    if not df.empty:
        x = df["mid_mb"].to_numpy()
        y = np.clip(df["value"].to_numpy(dtype=float), 0, ymax)
        if fill:
            ax.fill_between(x, 0, y, color=color, alpha=0.95, linewidth=0)
        else:
            ax.vlines(x, 0, y, color=color, lw=0.35)
    ax.set_ylim(0, ymax)
    ax.set_ylabel("")
    ax.text(0.012, 0.80, label, transform=ax.transAxes, fontsize=20, fontweight="bold", family="Arial")
    ax.tick_params(axis="y", labelsize=20, length=4)
    ax.tick_params(axis="x", labelbottom=False, length=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.9)
    ax.spines["bottom"].set_linewidth(0.9)


def main():
    args = parse_args()
    if args.font_file:
        font_path = Path(args.font_file)
        if not font_path.exists():
            raise SystemExit(f"ERROR: font file not found: {font_path}")
        font_manager.fontManager.addfont(str(font_path))
        font_name = font_manager.FontProperties(fname=str(font_path)).get_name()
        plt.rcParams.update({"font.family": font_name})

    try:
        import cooler
    except ImportError:
        raise SystemExit("ERROR: Python package 'cooler' is required. Activate the Hi-C conda environment first.")

    clr = cooler.Cooler(args.cool)
    chrom = args.chromosome
    if chrom not in clr.chromnames:
        raise SystemExit(f"ERROR: {chrom} not found in {args.cool}")
    chrom_len = int(clr.chromsizes[chrom])
    binsize = int(clr.binsize)
    mat = np.asarray(clr.matrix(balance=False).fetch(chrom), dtype=float)
    if args.hic_bin_step < 1:
        raise SystemExit("ERROR: --hic-bin-step must be >= 1")
    if args.hic_bin_step > 1:
        mat = mat[::args.hic_bin_step, ::args.hic_bin_step]
        binsize = binsize * args.hic_bin_step
    mat[~np.isfinite(mat)] = np.nan
    img = np.log1p(mat)
    finite = img[np.isfinite(img)]
    vmax = np.nanpercentile(finite, args.hic_vmax_percentile) if finite.size else 1
    if not np.isfinite(vmax) or vmax <= 0:
        vmax = 1

    n = img.shape[0]
    edges = np.arange(n + 1, dtype=float) * binsize
    edges[-1] = chrom_len
    x_edges = (edges[:, None] + edges[None, :]) / 2 / 1e6
    y_edges = (edges[None, :] - edges[:, None]) / 2 / 1e6
    tri_img = np.ma.array(img, mask=np.tril(np.ones_like(img, dtype=bool), k=-1))
    max_y = args.max_distance / 1e6 if args.max_distance else chrom_len / 2 / 1e6
    xlim = (0, chrom_len / 1e6)

    cen = read_centromere(args.centromere, chrom)
    copia, gypsy, satdna, gene, cg, chg, chh = scale_mtg_tracks(
        read_track(args.copia, chrom),
        read_track(args.gypsy, chrom),
        read_track(args.satdna, chrom),
        read_track(args.gene, chrom),
        read_track(args.cg, chrom),
        read_track(args.chg, chrom),
        read_track(args.chh, chrom),
    )
    tracks = [
        (copia, "#ff0000", "Copia", 1.0),
        (gypsy, "#00b7eb", "Gypsy", 1.0),
        (satdna, "#7b3294", "SatDNA", 1.0),
        (gene, "#16851c", "Genes", 1.0),
        (read_expr(args.expression, chrom, "leaf"), "#19a64a", "Leaf expression", 100.0),
        (read_expr(args.expression, chrom, "fruit"), "#f28c28", "Fruit expression", 100.0),
        (cg, "#d90000", "CG", 1.0),
        (chg, "#001fff", "CHG", 1.0),
        (chh, "#f79646", "CHH", 0.5),
    ]
    ltr = read_ticks(args.gypsy_ltr, chrom)
    intr = read_ticks(args.gypsy_int, chrom)

    heights = [3.2] + [0.82] * len(tracks) + [0.35, 0.35]
    fig = plt.figure(figsize=(8.0, 15.0), dpi=args.dpi)
    gs = fig.add_gridspec(len(heights), 1, height_ratios=heights, hspace=0.08)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.pcolormesh(x_edges, y_edges, tri_img, cmap=red_cmap(), vmin=0, vmax=vmax, shading="flat", rasterized=True)
    ax0.set_xlim(*xlim)
    ax0.set_ylim(0, max_y)
    ax0.set_aspect("equal")
    ax0.set_xticks([])
    ax0.set_yticks([])
    for sp in ax0.spines.values():
        sp.set_visible(False)

    axes = []
    for i, (df, color, label, ymax) in enumerate(tracks, start=1):
        ax = fig.add_subplot(gs[i, 0], sharex=ax0)
        draw_track(ax, df, color, label, ymax, cen=cen)
        ax.set_xlim(*xlim)
        axes.append(ax)

    for j, (tick_df, color) in enumerate([(ltr, "#001fff"), (intr, "#00b7eb")], start=1 + len(tracks)):
        ax = fig.add_subplot(gs[j, 0], sharex=ax0)
        if cen:
            ax.axvline(cen[0], color="black", lw=1.0, ls=(0, (4, 4)))
            ax.axvline(cen[1], color="black", lw=1.0, ls=(0, (4, 4)))
        if not tick_df.empty:
            ax.vlines(tick_df["mid_mb"].to_numpy(), 0, 1, color=color, lw=0.25, alpha=0.8)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
        ax.tick_params(axis="x", labelbottom=False, length=0)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        axes.append(ax)

    axes[-1].tick_params(axis="x", labelbottom=True, length=4, labelsize=20)
    axes[-1].set_xlabel(f"Position on {chrom} (Mb)", fontsize=20, family="Arial")
    fig.subplots_adjust(left=0.09, right=0.99, top=0.995, bottom=0.045)

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_prefix.with_suffix(".png"), dpi=args.dpi, facecolor="white")
    if not args.no_pdf:
        fig.savefig(out_prefix.with_suffix(".pdf"), facecolor="white")
    print(out_prefix.with_suffix(".png"))
    if not args.no_pdf:
        print(out_prefix.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
