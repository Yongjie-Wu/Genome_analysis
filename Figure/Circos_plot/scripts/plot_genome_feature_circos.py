#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def parse_args():
    parser = argparse.ArgumentParser(
        description="Draw a multi-track genome feature circos plot from genome-wide window track files."
    )
    parser.add_argument("--fai", required=True, help="Genome FASTA index file.")
    parser.add_argument("--centromere", required=True, help="Centromere range file: chr start end.")
    parser.add_argument("--gypsy", required=True, help="Genome-wide Gypsy density track: chr start end value.")
    parser.add_argument("--copia", required=True, help="Genome-wide Copia density track: chr start end value.")
    parser.add_argument("--gene", required=True, help="Genome-wide Gene density track: chr start end value.")
    parser.add_argument("--satdna", required=True, help="Genome-wide SatDNA density track: chr start end value.")
    parser.add_argument("--cg", required=True, help="Genome-wide CG methylation track: chr start end value.")
    parser.add_argument("--chg", required=True, help="Genome-wide CHG methylation track: chr start end value.")
    parser.add_argument("--chh", required=True, help="Genome-wide CHH methylation track: chr start end value.")
    parser.add_argument("--out-prefix", required=True, help="Output prefix without extension.")
    parser.add_argument("--expression", default=None, help="Optional expression file with chr/start/end/leaf/fruit.")
    parser.add_argument("--atac", default=None, help="Optional ATAC window TSV with chr/start/end/CPM_per_chr_mapped.")
    parser.add_argument("--chr-regex", default=r"^Chr[0-9]+A$", help="Chromosome names to plot.")
    parser.add_argument("--chr-suffix", default="A", help="Suffix removed from labels. Use none to keep labels.")
    parser.add_argument("--gap-bp", type=int, default=1500000, help="Pseudo gap size between chromosomes.")
    parser.add_argument("--plot-span-deg", type=float, default=350, help="Angular span in degrees.")
    parser.add_argument("--dpi", type=int, default=300, help="PNG resolution.")
    return parser.parse_args()


def read_chr_info(fai, chr_regex):
    chr_df = pd.read_csv(fai, sep="\t", header=None, usecols=[0, 1], names=["chr", "length"])
    rx = re.compile(chr_regex)
    chr_df = chr_df[chr_df["chr"].astype(str).map(lambda x: bool(rx.match(x)))].copy()
    if chr_df.empty:
        raise SystemExit(f"ERROR: no chromosomes matched --chr-regex {chr_regex}")

    extracted = chr_df["chr"].str.extract(r"(\d+)")[0]
    if extracted.notna().all():
        chr_df["ord"] = extracted.astype(int)
        chr_df = chr_df.sort_values("ord")
    else:
        chr_df = chr_df.sort_values("chr")
    return chr_df.reset_index(drop=True)


def prepare_offsets(chr_df, gap_bp):
    starts = []
    offset = 0
    for length in chr_df["length"]:
        starts.append(offset)
        offset += int(length) + gap_bp
    chr_df = chr_df.copy()
    chr_df["offset"] = starts
    return chr_df, offset


def read_track_file(path, chroms, value_col=None):
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: track file not found or empty: {path}")
    try:
        df = pd.read_csv(path, sep="\t")
        if {"chr", "start", "end"}.issubset(df.columns):
            if value_col is None:
                candidates = [c for c in df.columns if c not in {"chr", "start", "end"}]
                if not candidates:
                    raise ValueError("no value column")
                value_col = candidates[0]
            df = df[["chr", "start", "end", value_col]].copy()
            df.columns = ["chr", "start", "end", "value"]
        else:
            raise ValueError("no header")
    except Exception:
        df = pd.read_csv(path, sep=r"\s+", header=None, usecols=[0, 1, 2, 3],
                         names=["chr", "start", "end", "value"])
    df = df[df["chr"].isin(chroms)].copy()
    df["start"] = pd.to_numeric(df["start"], errors="coerce")
    df["end"] = pd.to_numeric(df["end"], errors="coerce")
    df["value"] = pd.to_numeric(df["value"], errors="coerce")
    return df.dropna(subset=["start", "end", "value"])[["chr", "start", "end", "value"]]


def read_expression(path, chroms, column):
    if not path:
        return pd.DataFrame(columns=["chr", "start", "end", "value"])
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: expression file not found: {path}")
    df = pd.read_csv(path, sep="\t")
    required = {"chr", "start", "end", column}
    if not required.issubset(df.columns):
        raise SystemExit(f"ERROR: expression file must contain columns: {', '.join(sorted(required))}")
    df = df[df["chr"].isin(chroms)][["chr", "start", "end", column]].copy()
    df["value"] = np.log10(pd.to_numeric(df[column], errors="coerce").fillna(0) + 1)
    return df.dropna(subset=["value"])[["chr", "start", "end", "value"]]


def read_atac(path, chroms):
    if not path:
        return pd.DataFrame(columns=["chr", "start", "end", "value"])
    path = Path(path)
    if not path.exists() or path.stat().st_size == 0:
        raise SystemExit(f"ERROR: ATAC file not found: {path}")
    df = pd.read_csv(path, sep="\t")
    required = {"chr", "start", "end", "CPM_per_chr_mapped"}
    if not required.issubset(df.columns):
        raise SystemExit(f"ERROR: ATAC file must contain columns: {', '.join(sorted(required))}")
    df = df[df["chr"].isin(chroms)][["chr", "start", "end", "CPM_per_chr_mapped"]].copy()
    df["value"] = pd.to_numeric(df["CPM_per_chr_mapped"], errors="coerce")
    return df.dropna(subset=["value"])[["chr", "start", "end", "value"]]


def main():
    args = parse_args()
    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    chr_df = read_chr_info(args.fai, args.chr_regex)
    if len(chr_df) != 17:
        print(f"WARNING: expected 17 chromosomes for this figure style, found {len(chr_df)}")
    chr_df, total = prepare_offsets(chr_df, args.gap_bp)
    chroms = chr_df["chr"].tolist()
    chr_offset = dict(zip(chr_df["chr"], chr_df["offset"]))

    plot_span = np.deg2rad(args.plot_span_deg)
    plot_gap = 2 * np.pi - plot_span
    plot_start = plot_gap / 2

    def theta_for(chrom, pos):
        frac = (chr_offset[chrom] + np.asarray(pos)) / total
        return plot_start + plot_span * frac

    tracks = [
        ("Copia", "Copia density", read_track_file(args.copia, chroms), "#ff4b4b", None),
        ("Gypsy", "Gypsy density", read_track_file(args.gypsy, chroms), "#00b7eb", None),
        ("SatDNA", "SatDNA density", read_track_file(args.satdna, chroms), "#8e44ad", None),
        ("ATAC", "ATAC signal CPM", read_atac(args.atac, chroms), "#6b6ecf", None),
        ("Gene", "Gene density", read_track_file(args.gene, chroms), "#1f9d3a", None),
        ("LeafExpr", "Leaf expression log10(RPKM+1)", read_expression(args.expression, chroms, "leaf"), "#f28e2b", None),
        ("FruitExpr", "Fruit expression log10(RPKM+1)", read_expression(args.expression, chroms, "fruit"), "#a65628", None),
        ("CG", "CG methylation", read_track_file(args.cg, chroms), "#e60000", (0, 1)),
        ("CHG", "CHG methylation", read_track_file(args.chg, chroms), "#0b3dff", (0, 1)),
        ("CHH", "CHH methylation", read_track_file(args.chh, chroms), "#f4a261", (0, 0.30)),
    ]
    tracks = [x for x in tracks if not x[2].empty]
    if not tracks:
        raise SystemExit("ERROR: no non-empty track files were found.")

    cent = pd.read_csv(args.centromere, sep=r"\s+", header=None, names=["chr", "start", "end"])
    cent = cent[cent["chr"].isin(chroms)].copy()

    fig = plt.figure(figsize=(12, 12), dpi=args.dpi)
    ax = fig.add_subplot(111, projection="polar")
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1.08)
    ax.set_axis_off()
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    outer_bottom = 0.965
    outer_height = 0.022
    palette = ["#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
               "#8491B4", "#91D1C2", "#DC0000", "#7E6148"]
    for i, row in chr_df.iterrows():
        chrom = row["chr"]
        center = theta_for(chrom, row["length"] / 2)
        width = plot_span * row["length"] / total
        ax.bar(center, outer_height, width=width, bottom=outer_bottom, align="center",
               color=palette[i % len(palette)], edgecolor="white", linewidth=0.35)
        label = chrom if args.chr_suffix == "none" else chrom.removesuffix(args.chr_suffix)
        ax.text(center, 1.028, label, fontsize=10, ha="center", va="center")

    for _, row in chr_df.iterrows():
        chrom = row["chr"]
        for pos in range(0, int(row["length"]) + 1, 5_000_000):
            th = theta_for(chrom, pos)
            ax.plot([th, th], [0.94, 0.958], color="#666666", lw=0.35, alpha=0.8)

    for _, row in cent.iterrows():
        chrom = row["chr"]
        center = theta_for(chrom, (row["start"] + row["end"]) / 2)
        width = plot_span * (row["end"] - row["start"]) / total
        ax.bar(center, outer_height * 1.08, width=width, bottom=outer_bottom - outer_height * 0.04,
               align="center", color="#111111", edgecolor="white", linewidth=0.25, zorder=5)

    track_height = 0.055
    gap = 0.012
    start_r = 0.85
    legend_handles = []
    roman_labels = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII"]

    for idx, (_, label, df, color, fixed_ylim) in enumerate(tracks):
        bottom = start_r - idx * (track_height + gap)
        top = bottom + track_height
        theta_grid = np.linspace(plot_start, plot_start + plot_span, 720)
        ax.plot(theta_grid, np.full(720, bottom), color="#EEEEEE", lw=0.45)
        ax.plot(theta_grid, np.full(720, top), color="#EEEEEE", lw=0.45)

        if fixed_ylim is None:
            ymax = float(df["value"].quantile(0.99))
            if not np.isfinite(ymax) or ymax <= 0:
                ymax = float(df["value"].max()) if len(df) else 1.0
            if not np.isfinite(ymax) or ymax <= 0:
                ymax = 1.0
            ymin = 0.0
        else:
            ymin, ymax = fixed_ylim

        for chrom in chroms:
            sub = df[df["chr"] == chrom].sort_values("start")
            if sub.empty:
                continue
            mid = (sub["start"].to_numpy(dtype=float) + sub["end"].to_numpy(dtype=float)) / 2
            val = sub["value"].to_numpy(dtype=float)
            val = np.clip((val - ymin) / (ymax - ymin), 0, 1)
            th = theta_for(chrom, mid)
            rr = bottom + val * track_height
            ax.fill_between(th, bottom, rr, color=color, alpha=0.85, linewidth=0)
            ax.plot(th, rr, color=color, lw=0.18, alpha=0.95)

        roman = roman_labels[idx] if idx < len(roman_labels) else str(idx + 1)
        ax.text(np.deg2rad(0), bottom + track_height / 2, roman, fontsize=8,
                fontweight="bold", ha="center", va="center", color="#222222")
        legend_handles.append(Patch(facecolor=color, edgecolor="none", label=f"{roman}  {label}"))

    ax.legend(handles=legend_handles, loc="center", frameon=False, fontsize=8, ncol=1, handlelength=1.1)
    for ext in ["png", "pdf", "svg"]:
        path = out_prefix.with_suffix(f".{ext}")
        fig.savefig(path, bbox_inches="tight", facecolor="white")
        print(path)


if __name__ == "__main__":
    main()
