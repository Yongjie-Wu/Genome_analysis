#!/usr/bin/env python3
import argparse
from pathlib import Path


def read_fai(genome):
    fai = Path(str(genome) + ".fai")
    if not fai.exists():
        raise SystemExit(f"ERROR: fasta index not found: {fai}. Run samtools faidx first.")
    chroms = []
    with fai.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                chroms.append((fields[0], int(fields[1])))
    return chroms


def init_bins(chroms, window):
    bins = {}
    order = []
    for chrom, length in chroms:
        for start in range(0, length, window):
            end = min(start + window, length)
            key = (chrom, start)
            bins[key] = [end, 0.0, 0.0]
            order.append(key)
    return bins, order


def add_cov(infile, bins, window):
    with open(infile) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[0]
            try:
                pos = int(float(fields[1]))
            except ValueError:
                continue
            start = (pos // window) * window
            key = (chrom, start)
            if key not in bins:
                continue

            if len(fields) >= 11:
                try:
                    cov = float(fields[9])
                    meth = float(fields[10])
                except ValueError:
                    continue
                if cov <= 0:
                    continue
                if meth > 1:
                    meth = meth / 100.0
                bins[key][1] += cov * meth
                bins[key][2] += cov
            else:
                try:
                    meth = float(fields[3])
                except ValueError:
                    continue
                if meth > 1:
                    meth = meth / 100.0
                bins[key][1] += meth
                bins[key][2] += 1.0


def write_bins(outfile, bins, order):
    with open(outfile, "w") as out:
        for chrom, start in order:
            end, meth_sum, cov_sum = bins[(chrom, start)]
            value = meth_sum / cov_sum if cov_sum > 0 else 0.0
            out.write(f"{chrom}\t{start}\t{end}\t{value:.6f}\t{cov_sum:.2f}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--cpg", required=True)
    parser.add_argument("--chg", required=True)
    parser.add_argument("--chh", required=True)
    parser.add_argument("--window", required=True, type=int)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    chroms = read_fai(args.genome)

    for label, infile in [("CG", args.cpg), ("CHG", args.chg), ("CHH", args.chh)]:
        bins, order = init_bins(chroms, args.window)
        add_cov(infile, bins, args.window)
        outfile = outdir / f"{args.sample}.{label}.{args.window}.bed"
        write_bins(outfile, bins, order)
        print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
