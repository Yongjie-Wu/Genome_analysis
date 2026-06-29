#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Prepare genome-wide window tracks for genome feature circos plot.")
    p.add_argument("--sample", required=True)
    p.add_argument("--fai", required=True)
    p.add_argument("--centromere", required=True)
    p.add_argument("--gypsy", required=True)
    p.add_argument("--copia", required=True)
    p.add_argument("--gene", required=True)
    p.add_argument("--satdna", required=True)
    p.add_argument("--cpg", required=True)
    p.add_argument("--chg", required=True)
    p.add_argument("--chh", required=True)
    p.add_argument("--rpkm", default=None)
    p.add_argument("--outdir", required=True)
    p.add_argument("--window", type=int, default=30000)
    p.add_argument("--chr-regex", default=r"^Chr[0-9]+A$")
    return p.parse_args()


def read_chr_info(fai, chr_regex):
    df = pd.read_csv(fai, sep="\t", header=None, usecols=[0, 1], names=["chr", "length"])
    rx = re.compile(chr_regex)
    df = df[df["chr"].astype(str).map(lambda x: bool(rx.match(x)))].copy()
    if df.empty:
        raise SystemExit(f"ERROR: no chromosomes matched {chr_regex}")
    num = df["chr"].str.extract(r"(\d+)")[0]
    if num.notna().all():
        df["ord"] = num.astype(int)
        df = df.sort_values("ord")
    else:
        df = df.sort_values("chr")
    return df[["chr", "length"]].reset_index(drop=True)


def make_windows(chr_info, window):
    rows = []
    for chrom, length in chr_info[["chr", "length"]].itertuples(index=False):
        for start in range(0, int(length), window):
            end = min(start + window, int(length))
            rows.append((chrom, start, end, end - start))
    return pd.DataFrame(rows, columns=["chr", "start", "end", "window_size"])


def read_bed(path, chroms):
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 3:
                f = line.rstrip("\n").split()
            if len(f) < 3 or f[0] not in chroms:
                continue
            try:
                rows.append((f[0], int(float(f[1])), int(float(f[2]))))
            except ValueError:
                continue
    return pd.DataFrame(rows, columns=["chr", "start", "end"])


def read_gene(path, chroms):
    path = Path(path)
    rows = []
    if path.suffix.lower() in {".gff", ".gff3"}:
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[0] not in chroms or f[2] != "gene":
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
                if len(f) < 3 or f[0] not in chroms:
                    continue
                gid = f[3] if len(f) >= 4 else f"gene_{len(rows)+1}"
                rows.append((f[0], int(float(f[1])), int(float(f[2])), gid))
    return pd.DataFrame(rows, columns=["chr", "start", "end", "gene_id"])


def interval_density(windows, bed_df):
    values = np.zeros(len(windows), dtype=float)
    if bed_df.empty:
        return values
    for chrom, widx in windows.groupby("chr").groups.items():
        wsub = windows.loc[widx]
        bsub = bed_df[bed_df["chr"] == chrom]
        if bsub.empty:
            continue
        for b in bsub.itertuples(index=False):
            hits = wsub[(wsub["end"] > b.start) & (wsub["start"] < b.end)]
            if hits.empty:
                continue
            ov = np.minimum(hits["end"].to_numpy(), b.end) - np.maximum(hits["start"].to_numpy(), b.start)
            values[hits.index.to_numpy()] += np.maximum(ov, 0)
    return values / windows["window_size"].to_numpy()


def gene_count_track(windows, gene_df):
    values = np.zeros(len(windows), dtype=int)
    if gene_df.empty:
        return values
    for chrom, widx in windows.groupby("chr").groups.items():
        wsub = windows.loc[widx]
        gsub = gene_df[gene_df["chr"] == chrom]
        for g in gsub.itertuples(index=False):
            hits = wsub[(wsub["end"] > g.start) & (wsub["start"] < g.end)]
            values[hits.index.to_numpy()] += 1
    return values


def methylation_track(windows, cov_file, chroms):
    meth_sum = np.zeros(len(windows), dtype=float)
    cov_sum = np.zeros(len(windows), dtype=float)
    index = {}
    for chrom, idxs in windows.groupby("chr").groups.items():
        starts = windows.loc[idxs, "start"].to_numpy()
        index[chrom] = (starts, idxs.to_numpy())
    with open(cov_file) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                f = line.rstrip("\n").split()
            if len(f) < 4 or f[0] not in chroms or f[0] not in index:
                continue
            try:
                pos = int(float(f[1]))
            except ValueError:
                continue
            starts, idxs = index[f[0]]
            rel = np.searchsorted(starts, pos, side="right") - 1
            if rel < 0 or rel >= len(idxs):
                continue
            out_i = idxs[rel]
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
                meth_sum[out_i] += cov * val
                cov_sum[out_i] += cov
            else:
                try:
                    val = float(f[3])
                except ValueError:
                    continue
                if val > 1:
                    val /= 100.0
                meth_sum[out_i] += val
                cov_sum[out_i] += 1
    return np.divide(meth_sum, cov_sum, out=np.zeros_like(meth_sum), where=cov_sum > 0)


def write_track(windows, values, outfile):
    out = windows[["chr", "start", "end"]].copy()
    out["value"] = values
    out.to_csv(outfile, sep="\t", index=False, header=False)


def expression_track(windows, gene_df, rpkm_file, outfile):
    expr = pd.read_csv(rpkm_file, sep="\t")
    required = {"ID", "leaf", "fruit"}
    if not required.issubset(expr.columns):
        raise SystemExit("ERROR: RPKM file must contain columns: ID, leaf, fruit")
    g = gene_df.merge(expr, left_on="gene_id", right_on="ID", how="left")
    g[["leaf", "fruit"]] = g[["leaf", "fruit"]].fillna(0)
    leaf_sum = np.zeros(len(windows), dtype=float)
    fruit_sum = np.zeros(len(windows), dtype=float)
    count = np.zeros(len(windows), dtype=float)
    for chrom, widx in windows.groupby("chr").groups.items():
        wsub = windows.loc[widx]
        gsub = g[g["chr"] == chrom]
        for row in gsub.itertuples(index=False):
            hits = wsub[(wsub["end"] > row.start) & (wsub["start"] < row.end)]
            idx = hits.index.to_numpy()
            leaf_sum[idx] += float(row.leaf)
            fruit_sum[idx] += float(row.fruit)
            count[idx] += 1
    out = windows[["chr", "start", "end"]].copy()
    out["leaf"] = np.divide(leaf_sum, count, out=np.zeros_like(leaf_sum), where=count > 0)
    out["fruit"] = np.divide(fruit_sum, count, out=np.zeros_like(fruit_sum), where=count > 0)
    out.to_csv(outfile, sep="\t", index=False)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    chr_info = read_chr_info(args.fai, args.chr_regex)
    chroms = set(chr_info["chr"])
    windows = make_windows(chr_info, args.window)
    windows[["chr", "start", "end"]].to_csv(outdir / f"{args.sample}.windows.{args.window}bp.bed", sep="\t", index=False, header=False)

    gypsy = read_bed(args.gypsy, chroms)
    copia = read_bed(args.copia, chroms)
    satdna = read_bed(args.satdna, chroms)
    gene = read_gene(args.gene, chroms)

    write_track(windows, interval_density(windows, gypsy), outdir / f"{args.sample}.Gypsy.{args.window}bp.txt")
    write_track(windows, interval_density(windows, copia), outdir / f"{args.sample}.Copia.{args.window}bp.txt")
    write_track(windows, gene_count_track(windows, gene), outdir / f"{args.sample}.Gene.{args.window}bp.txt")
    write_track(windows, interval_density(windows, satdna), outdir / f"{args.sample}.SatDNA.{args.window}bp.txt")
    write_track(windows, methylation_track(windows, args.cpg, chroms), outdir / f"{args.sample}.CG.{args.window}bp.txt")
    write_track(windows, methylation_track(windows, args.chg, chroms), outdir / f"{args.sample}.CHG.{args.window}bp.txt")
    write_track(windows, methylation_track(windows, args.chh, chroms), outdir / f"{args.sample}.CHH.{args.window}bp.txt")

    if args.rpkm:
        expression_track(windows, gene, args.rpkm, outdir / f"{args.sample}.leaf_fruit_expression.{args.window}bp.tsv")

    print(f"Wrote circos inputs to {outdir}")


if __name__ == "__main__":
    main()
