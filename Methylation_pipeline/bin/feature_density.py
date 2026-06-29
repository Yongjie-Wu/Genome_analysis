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


def gene_bed_from_gff3(gff3, outfile):
    kept = 0
    with open(gff3) as inp, open(outfile, "w") as out:
        for line in inp:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2].lower()
            if feature not in {"gene", "mrna", "transcript"}:
                continue
            try:
                start = max(0, int(fields[3]) - 1)
                end = int(fields[4])
            except ValueError:
                continue
            out.write(f"{fields[0]}\t{start}\t{end}\t{fields[8]}\t0\t{fields[6]}\n")
            kept += 1
    if kept == 0:
        raise SystemExit(f"ERROR: no gene/mRNA/transcript records parsed from {gff3}")


def init_windows(chroms, window):
    rows = []
    index = {}
    for chrom, length in chroms:
        for start in range(0, length, window):
            end = min(start + window, length)
            key = (chrom, start)
            index[key] = len(rows)
            rows.append({
                "chr": chrom,
                "start": start,
                "end": end,
                "gene_density": 0,
                "satdna_density": 0,
                "copia_density": 0,
                "gypsy_density": 0,
                "in_centromere": 0,
            })
    return rows, index


def add_bed_density(path, rows, index, window, column):
    if path == "none":
        return
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return
    with p.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            try:
                start = int(float(fields[1]))
                end = int(float(fields[2]))
            except ValueError:
                continue
            if end <= start:
                continue
            first = (start // window) * window
            last = ((end - 1) // window) * window
            for win_start in range(first, last + 1, window):
                key = (chrom, win_start)
                if key not in index:
                    continue
                row = rows[index[key]]
                ov_start = max(start, row["start"])
                ov_end = min(end, row["end"])
                if ov_end > ov_start:
                    row[column] += ov_end - ov_start


def add_centromere(path, rows, index, window):
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split()
            if len(fields) < 3:
                continue
            chrom = fields[0]
            try:
                start = int(float(fields[1]))
                end = int(float(fields[2]))
            except ValueError:
                continue
            first = (start // window) * window
            last = ((end - 1) // window) * window
            for win_start in range(first, last + 1, window):
                key = (chrom, win_start)
                if key not in index:
                    continue
                row = rows[index[key]]
                ov_start = max(start, row["start"])
                ov_end = min(end, row["end"])
                if ov_end > ov_start:
                    row["in_centromere"] = 1


def normalize(rows):
    for row in rows:
        length = row["end"] - row["start"]
        if length <= 0:
            continue
        for col in ["gene_density", "satdna_density", "copia_density", "gypsy_density"]:
            row[col] = row[col] / length


def write_table(rows, outfile):
    cols = ["chr", "start", "end", "gene_density", "satdna_density", "copia_density", "gypsy_density", "in_centromere"]
    with open(outfile, "w") as out:
        out.write("\t".join(cols) + "\n")
        for row in rows:
            out.write("\t".join(str(row[c]) for c in cols) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--gff3", required=True)
    parser.add_argument("--centromere", required=True)
    parser.add_argument("--satdna", required=True)
    parser.add_argument("--copia", required=True)
    parser.add_argument("--gypsy", required=True)
    parser.add_argument("--window", type=int, required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gene_bed = outdir / f"{args.sample}.gene.bed"
    gene_bed_from_gff3(args.gff3, gene_bed)

    rows, index = init_windows(read_fai(args.genome), args.window)
    add_bed_density(str(gene_bed), rows, index, args.window, "gene_density")
    add_bed_density(args.satdna, rows, index, args.window, "satdna_density")
    add_bed_density(args.copia, rows, index, args.window, "copia_density")
    add_bed_density(args.gypsy, rows, index, args.window, "gypsy_density")
    add_centromere(args.centromere, rows, index, args.window)
    normalize(rows)

    outfile = outdir / f"{args.sample}.feature_density.{args.window}.tsv"
    write_table(rows, outfile)
    print(f"Wrote {gene_bed}")
    print(f"Wrote {outfile}")


if __name__ == "__main__":
    main()
