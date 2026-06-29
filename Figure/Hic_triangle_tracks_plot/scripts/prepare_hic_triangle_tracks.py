#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Prepare one-chromosome tracks for Hi-C triangle plus multi-track plot.")
    p.add_argument("--sample", required=True)
    p.add_argument("--chromosome", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--centromere", required=True)
    p.add_argument("--gypsy", required=True)
    p.add_argument("--copia", required=True)
    p.add_argument("--gene", required=True)
    p.add_argument("--satdna", required=True)
    p.add_argument("--cpg", required=True)
    p.add_argument("--chg", required=True)
    p.add_argument("--chh", required=True)
    p.add_argument("--rpkm")
    p.add_argument("--gypsy-ltr")
    p.add_argument("--gypsy-int")
    p.add_argument("--outdir", required=True)
    p.add_argument("--window", type=int, default=30000)
    return p.parse_args()


def chrom_length(fai, chrom):
    with open(fai) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2 and f[0] == chrom:
                return int(f[1])
    raise SystemExit(f"ERROR: chromosome {chrom} not found in {fai}")


def make_windows(chrom, length, window):
    rows = []
    for start in range(0, length, window):
        end = min(start + window, length)
        rows.append((chrom, start, end, end - start))
    return pd.DataFrame(rows, columns=["chr", "start", "end", "window_size"])


def read_bed(path, chrom):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                f = line.rstrip("\n").split()
            if len(f) < 3 or f[0] != chrom:
                continue
            try:
                rows.append((f[0], int(float(f[1])), int(float(f[2]))))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["chr", "start", "end"])


def read_gene(path, chrom):
    path = Path(path)
    rows = []
    if path.suffix.lower() in {".gff", ".gff3"}:
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[0] != chrom or f[2] != "gene":
                    continue
                m = re.search(r"ID=([^;]+)", f[8])
                gid = m.group(1) if m else f"gene_{len(rows)+1}"
                rows.append((f[0], int(f[3]) - 1, int(f[4]), gid))
    else:
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 3:
                    f = line.rstrip("\n").split()
                if len(f) < 3 or f[0] != chrom:
                    continue
                gid = f[3] if len(f) >= 4 else f"gene_{len(rows)+1}"
                rows.append((f[0], int(float(f[1])), int(float(f[2])), gid))
    return pd.DataFrame(rows, columns=["chr", "start", "end", "gene_id"])


def interval_density(windows, bed_df):
    values = np.zeros(len(windows), dtype=float)
    if bed_df.empty:
        return values
    for b in bed_df.itertuples(index=False):
        hits = windows[(windows["end"] > b.start) & (windows["start"] < b.end)]
        if hits.empty:
            continue
        ov = np.minimum(hits["end"].to_numpy(), b.end) - np.maximum(hits["start"].to_numpy(), b.start)
        values[hits.index.to_numpy()] += np.maximum(ov, 0)
    return values / windows["window_size"].to_numpy()


def gene_count(windows, gene_df):
    values = np.zeros(len(windows), dtype=int)
    for g in gene_df.itertuples(index=False):
        hits = windows[(windows["end"] > g.start) & (windows["start"] < g.end)]
        values[hits.index.to_numpy()] += 1
    return values


def methylation_track(windows, cov_file, chrom):
    meth_sum = np.zeros(len(windows), dtype=float)
    cov_sum = np.zeros(len(windows), dtype=float)
    starts = windows["start"].to_numpy()
    with open(cov_file) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                f = line.rstrip("\n").split()
            if len(f) < 4 or f[0] != chrom:
                continue
            try:
                pos = int(float(f[1]))
            except ValueError:
                continue
            idx = np.searchsorted(starts, pos, side="right") - 1
            if idx < 0 or idx >= len(windows):
                continue
            if len(f) >= 11:
                try:
                    cov = float(f[9])
                    val = float(f[10])
                except ValueError:
                    continue
                if cov <= 0:
                    continue
                if val > 1:
                    val /= 100.0
                meth_sum[idx] += cov * val
                cov_sum[idx] += cov
            else:
                try:
                    val = float(f[3])
                except ValueError:
                    continue
                if val > 1:
                    val /= 100.0
                meth_sum[idx] += val
                cov_sum[idx] += 1
    return np.divide(meth_sum, cov_sum, out=np.zeros_like(meth_sum), where=cov_sum > 0)


def write_track(windows, values, path):
    out = windows[["chr", "start", "end"]].copy()
    out["value"] = values
    out.to_csv(path, sep="\t", index=False, header=False)


def expression_track(windows, gene_df, rpkm_file, path):
    expr = pd.read_csv(rpkm_file, sep="\t")
    required = {"ID", "leaf", "fruit"}
    if not required.issubset(expr.columns):
        raise SystemExit("ERROR: RPKM file must contain columns: ID, leaf, fruit")
    g = gene_df.merge(expr, left_on="gene_id", right_on="ID", how="left")
    g[["leaf", "fruit"]] = g[["leaf", "fruit"]].fillna(0)
    leaf_sum = np.zeros(len(windows), dtype=float)
    fruit_sum = np.zeros(len(windows), dtype=float)
    count = np.zeros(len(windows), dtype=float)
    for row in g.itertuples(index=False):
        hits = windows[(windows["end"] > row.start) & (windows["start"] < row.end)]
        idx = hits.index.to_numpy()
        leaf_sum[idx] += float(row.leaf)
        fruit_sum[idx] += float(row.fruit)
        count[idx] += 1
    out = windows[["chr", "start", "end"]].copy()
    out["leaf"] = np.divide(leaf_sum, count, out=np.zeros_like(leaf_sum), where=count > 0)
    out["fruit"] = np.divide(fruit_sum, count, out=np.zeros_like(fruit_sum), where=count > 0)
    out.to_csv(path, sep="\t", index=False)


def copy_chr_bed(infile, chrom, outfile):
    df = read_bed(infile, chrom)
    if df.empty:
        Path(outfile).write_text("")
    else:
        df.to_csv(outfile, sep="\t", index=False, header=False)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    length = chrom_length(args.fai, args.chromosome)
    windows = make_windows(args.chromosome, length, args.window)
    prefix = f"{args.sample}.{args.chromosome}"
    windows[["chr", "start", "end"]].to_csv(outdir / f"{prefix}.windows.{args.window}bp.bed", sep="\t", index=False, header=False)

    gypsy = read_bed(args.gypsy, args.chromosome)
    copia = read_bed(args.copia, args.chromosome)
    satdna = read_bed(args.satdna, args.chromosome)
    gene = read_gene(args.gene, args.chromosome)

    write_track(windows, interval_density(windows, copia), outdir / f"{prefix}.Copia.{args.window}bp.txt")
    write_track(windows, interval_density(windows, gypsy), outdir / f"{prefix}.Gypsy.{args.window}bp.txt")
    write_track(windows, interval_density(windows, satdna), outdir / f"{prefix}.SatDNA.{args.window}bp.txt")
    write_track(windows, interval_density(windows, gene[["chr", "start", "end"]]), outdir / f"{prefix}.Gene.{args.window}bp.txt")
    write_track(windows, methylation_track(windows, args.cpg, args.chromosome), outdir / f"{prefix}.CG.{args.window}bp.txt")
    write_track(windows, methylation_track(windows, args.chg, args.chromosome), outdir / f"{prefix}.CHG.{args.window}bp.txt")
    write_track(windows, methylation_track(windows, args.chh, args.chromosome), outdir / f"{prefix}.CHH.{args.window}bp.txt")

    if args.rpkm:
        expression_track(windows, gene, args.rpkm, outdir / f"{prefix}.leaf_fruit_expression.{args.window}bp.tsv")
    if args.gypsy_ltr:
        copy_chr_bed(args.gypsy_ltr, args.chromosome, outdir / f"{prefix}.Gypsy_LTR.bed")
    if args.gypsy_int:
        copy_chr_bed(args.gypsy_int, args.chromosome, outdir / f"{prefix}.Gypsy_INT.bed")

    print(f"Wrote tracks to {outdir}")


if __name__ == "__main__":
    main()
