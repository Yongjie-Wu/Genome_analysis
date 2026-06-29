#!/usr/bin/env python3
"""
Usage: python bin/atac_window_features.py <config.sh> [hap1|hap2|both]

Creates 09.atac_windows/<hap>/ tables and plots:
- <project>.<hap>.allchr.<window>bp.ATAC_SatDNA_gene_expression.tsv
- <project>.<hap>.centromere_vs_arm.summary.tsv
- <project>.<hap>.ATAC_correlations.tsv
- <project>.<hap>.ATAC_correlations_scatter.png/pdf
- <project>.<hap>.centromere_vs_arm.*.png/pdf
"""
from __future__ import annotations
import os
import re
import sys
import subprocess
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def parse_config(path: Path) -> dict[str, str]:
    cfg = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        k = k.strip()
        v = v.strip().strip('"').strip("'")
        cfg[k] = v
    return cfg


def run(cmd):
    print("[cmd]", " ".join(map(str, cmd)), flush=True)
    subprocess.run(list(map(str, cmd)), check=True)


def spearman(x, y):
    tmp = pd.DataFrame({"x": x, "y": y}).replace([np.inf, -np.inf], np.nan).dropna()
    if len(tmp) < 3 or tmp.x.nunique() < 2 or tmp.y.nunique() < 2:
        return np.nan
    return tmp.x.rank().corr(tmp.y.rank())


def read_expression(path: Path) -> pd.DataFrame:
    expr = pd.read_csv(path, sep="\t")
    if "ID" not in expr.columns:
        expr = expr.rename(columns={expr.columns[0]: "ID"})
    numeric_cols = [c for c in expr.columns if c != "ID"]
    for c in numeric_cols:
        expr[c] = pd.to_numeric(expr[c], errors="coerce").fillna(0)
    if "leaf" not in expr.columns or "fruit" not in expr.columns:
        expr["mean_RPKM"] = expr[numeric_cols].mean(axis=1) if numeric_cols else 0
        expr["leaf"] = expr["mean_RPKM"]
        expr["fruit"] = expr["mean_RPKM"]
    else:
        expr["mean_RPKM"] = (expr["leaf"] + expr["fruit"]) / 2
    return expr[["ID", "leaf", "fruit", "mean_RPKM"]]


def analyze_hap(cfg: dict[str, str], hap: str):
    project = cfg["PROJECT_NAME"]
    project_dir = Path(cfg["PROJECT_DIR"])
    window = int(cfg.get("WINDOW_SIZE", "50000"))
    fa = Path(cfg[f"HAP{hap[-1]}_FA"])
    gff = Path(cfg[f"HAP{hap[-1]}_GFF"])
    cent = Path(cfg.get(f"HAP{hap[-1]}_CENTROMERE", ""))
    sat_path = Path(cfg.get(f"HAP{hap[-1]}_SATDNA_BED", ""))
    rpkm_path = Path(cfg.get(f"HAP{hap[-1]}_RPKM", ""))
    bam = project_dir / f"03.bam_filter/{project}.{hap}.merge.rmdup.bam"

    required = [fa, Path(str(fa) + ".fai"), gff, cent, sat_path, rpkm_path, bam]
    missing = [str(p) for p in required if not str(p) or not p.exists() or p.stat().st_size == 0]
    if missing:
        print(f"[skip] {hap}: missing optional step09 input(s):")
        for m in missing:
            print("  ", m)
        return

    out = project_dir / "09.atac_windows" / hap
    feature = out / "features"
    plots = out / "plots"
    tmp = out / "tmp"
    win_dir = out / "windows"
    for d in [feature, plots, tmp, win_dir]:
        d.mkdir(parents=True, exist_ok=True)

    chr_info = pd.read_csv(fa.with_suffix(fa.suffix + ".fai"), sep="\t", header=None, usecols=[0, 1], names=["chr", "length"])
    chr_info = chr_info[chr_info["chr"].str.match(r"^Chr\d+[A-Za-z]?$", na=False)].copy()
    if chr_info.empty:
        chr_info = pd.read_csv(fa.with_suffix(fa.suffix + ".fai"), sep="\t", header=None, usecols=[0, 1], names=["chr", "length"])
    chr_info["ord"] = range(len(chr_info))
    chroms = chr_info["chr"].tolist()

    cent_df = pd.read_csv(cent, sep=r"\s+", header=None, names=["chr", "cen_start", "cen_end"])
    cent_map = {r.chr: (int(r.cen_start), int(r.cen_end)) for r in cent_df.itertuples(index=False)}
    sat = pd.read_csv(sat_path, sep="\t", header=None, usecols=[0, 1, 2], names=["chr", "start", "end"])
    expr = read_expression(rpkm_path)

    genes = []
    with open(gff) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            m = re.search(r"ID=([^;]+)", parts[8])
            gid = m.group(1) if m else f"gene_{len(genes)+1}"
            genes.append((parts[0], int(parts[3]) - 1, int(parts[4]), gid, parts[6]))
    gene = pd.DataFrame(genes, columns=["chr", "start", "end", "gene_id", "strand"])
    gene = gene.merge(expr, left_on="gene_id", right_on="ID", how="left")
    gene[["leaf", "fruit", "mean_RPKM"]] = gene[["leaf", "fruit", "mean_RPKM"]].fillna(0)
    gene.to_csv(feature / f"{project}.{hap}.genes.with_expression.tsv", sep="\t", index=False)

    all_rows = []
    for row in chr_info.itertuples(index=False):
        chr_name = row.chr
        length = int(row.length)
        print(f"[info] {hap} {chr_name}", flush=True)
        win = win_dir / f"{hap}.{chr_name}.{window}bp.windows.bed"
        with open(win, "w") as outfh:
            for s in range(0, length, window):
                e = min(s + window, length)
                outfh.write(f"{chr_name}\t{s}\t{e}\n")
        chr_bam = tmp / f"{project}.{hap}.{chr_name}.merge.rmdup.bam"
        if not chr_bam.exists() or chr_bam.stat().st_size == 0 or not Path(str(chr_bam) + ".bai").exists():
            run(["samtools", "view", "-@", "4", "-b", bam, chr_name, "-o", chr_bam])
            run(["samtools", "index", "-@", "4", chr_bam])
        mapped_chr = int(subprocess.check_output(["samtools", "view", "-@", "4", "-c", "-F", "4", str(chr_bam)]).decode().strip())
        p = subprocess.run(["bedtools", "coverage", "-a", str(win), "-b", str(chr_bam), "-counts"], check=True, capture_output=True, text=True)
        rows = []
        for line in p.stdout.splitlines():
            c, s, e, count = line.split("\t")[:4]
            s, e, count = int(s), int(e), int(count)
            w = e - s
            rows.append((c, s, e, w, count, count / (w / 1000), count * 1_000_000 / mapped_chr if mapped_chr else np.nan))
        df = pd.DataFrame(rows, columns=["chr", "start", "end", "window_size", "read_count", "reads_per_kb", "CPM_per_chr_mapped"])
        cs, ce = cent_map.get(chr_name, (-1, -1))
        df["region"] = np.where((df["end"] > cs) & (df["start"] < ce), "centromere", "arm")
        df["mid"] = (df["start"] + df["end"]) / 2
        ss = sat[sat["chr"] == chr_name]
        df["SatDNA_bp"] = 0
        df["SatDNA_count"] = 0
        for i, w in df.iterrows():
            ov = ss[(ss.end > w.start) & (ss.start < w.end)]
            if not ov.empty:
                overlap = np.minimum(ov.end.to_numpy(), w.end) - np.maximum(ov.start.to_numpy(), w.start)
                overlap = overlap[overlap > 0]
                df.at[i, "SatDNA_bp"] = int(overlap.sum())
                df.at[i, "SatDNA_count"] = len(overlap)
        df["SatDNA_fraction"] = df["SatDNA_bp"] / df["window_size"]
        gg_chr = gene[gene["chr"] == chr_name]
        df["gene_count"] = 0
        df["expressed_gene_count"] = 0
        df["leaf_mean_RPKM"] = 0.0
        df["fruit_mean_RPKM"] = 0.0
        df["mean_RPKM"] = 0.0
        for i, w in df.iterrows():
            gg = gg_chr[(gg_chr.end > w.start) & (gg_chr.start < w.end)]
            df.at[i, "gene_count"] = len(gg)
            if len(gg):
                df.at[i, "leaf_mean_RPKM"] = gg.leaf.mean()
                df.at[i, "fruit_mean_RPKM"] = gg.fruit.mean()
                df.at[i, "mean_RPKM"] = gg.mean_RPKM.mean()
                df.at[i, "expressed_gene_count"] = int((gg.mean_RPKM > 1).sum())
        df["log10_mean_RPKM_plus1"] = np.log10(df["mean_RPKM"] + 1)
        all_rows.append(df)

    all_df = pd.concat(all_rows, ignore_index=True)
    all_table = feature / f"{project}.{hap}.allchr.{window}bp.ATAC_SatDNA_gene_expression.tsv"
    all_df.to_csv(all_table, sep="\t", index=False)

    summary = []
    for (chr_name, region), dd in all_df.groupby(["chr", "region"]):
        summary.append({"chr": chr_name, "region": region, "n_windows": len(dd), "ATAC_mean": dd.read_count.mean(), "ATAC_median": dd.read_count.median(), "SatDNA_fraction_mean": dd.SatDNA_fraction.mean(), "gene_count_mean": dd.gene_count.mean(), "mean_RPKM_mean": dd.mean_RPKM.mean()})
    for region, dd in all_df.groupby("region"):
        summary.append({"chr": "ALL", "region": region, "n_windows": len(dd), "ATAC_mean": dd.read_count.mean(), "ATAC_median": dd.read_count.median(), "SatDNA_fraction_mean": dd.SatDNA_fraction.mean(), "gene_count_mean": dd.gene_count.mean(), "mean_RPKM_mean": dd.mean_RPKM.mean()})
    pd.DataFrame(summary).to_csv(feature / f"{project}.{hap}.centromere_vs_arm.summary.tsv", sep="\t", index=False)

    corr_rows = []
    for scope, dd in [("ALL", all_df)] + [(r, d) for r, d in all_df.groupby("region")]:
        for ycol in ["SatDNA_fraction", "gene_count", "mean_RPKM", "leaf_mean_RPKM", "fruit_mean_RPKM", "log10_mean_RPKM_plus1"]:
            corr_rows.append({"scope": scope, "x": "ATAC_read_count", "y": ycol, "spearman_rho": spearman(dd.read_count, dd[ycol]), "n_windows": len(dd)})
    pd.DataFrame(corr_rows).to_csv(feature / f"{project}.{hap}.ATAC_correlations.tsv", sep="\t", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4), dpi=300)
    specs = [("SatDNA_fraction", "SatDNA fraction", "#8e44ad"), ("log10_mean_RPKM_plus1", "Expression log10(RPKM+1)", "#f28e2b"), ("gene_count", "Gene count / window", "#1f9d3a")]
    for ax, (ycol, ylabel, color) in zip(axes, specs):
        for region, marker, c in [("arm", "o", color), ("centromere", "s", "#111111")]:
            dd = all_df[all_df.region == region]
            ax.scatter(dd[ycol], dd.read_count, s=8, alpha=0.35, marker=marker, color=c, edgecolor="none", label=region)
        rho = spearman(all_df.read_count, all_df[ycol])
        ax.set_title(f"rho={rho:.3f}", fontsize=10)
        ax.set_xlabel(ylabel)
        ax.set_ylabel(f"ATAC reads / {window} bp")
        ax.grid(color="#dddddd", linewidth=0.5, alpha=0.8)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(plots / f"{project}.{hap}.allchr.{window}bp.ATAC_correlations_scatter.png", bbox_inches="tight", facecolor="white")
    fig.savefig(plots / f"{project}.{hap}.allchr.{window}bp.ATAC_correlations_scatter.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("[done]", all_table)


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    cfg = parse_config(Path(sys.argv[1]))
    mode = sys.argv[2] if len(sys.argv) > 2 else "both"
    if mode not in {"hap1", "hap2", "both"}:
        raise SystemExit("mode must be hap1, hap2, or both")
    haps = ["hap1", "hap2"] if mode == "both" else [mode]
    for hap in haps:
        analyze_hap(cfg, hap)

if __name__ == "__main__":
    main()
