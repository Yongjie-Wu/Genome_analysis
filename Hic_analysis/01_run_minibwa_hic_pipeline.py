#!/usr/bin/env python3
"""
Run the minibwa Hi-C pipeline to generate .cool and .hic files.

Usage:
    python run_minibwa_hic_pipeline.py <hic_file1> <hic_file2> <genome.fasta>

Example:
    python run_minibwa_hic_pipeline.py B186duli-ye_L8_310X10.R1.fastq.gz B186duli-ye_L8_310X10.R2.fastq.gz Pyrbe_duli_hap2.fa

Outputs are written to the current working directory. The output prefix is
derived from the genome FASTA name:
    <genome_basename>.minibwa.10k.cool
    <genome_basename>.minibwa.hic
"""

from __future__ import annotations

import os
import shlex
import subprocess
import sys
from pathlib import Path


THREADS = int(os.environ.get("THREADS", "40"))
MINIBWA = Path(os.environ.get("MINIBWA", "/storage3/wuyongjie/software/minibwa/minibwa"))
JUICER_TOOLS = Path(os.environ.get("JUICER_TOOLS", "/storage3/wuyongjie/software/juicer_tools.jar"))
HIC_AB_BIN = Path(os.environ.get("HIC_AB_BIN", "/storage3/wuyongjie/miniconda3/envs/hic_ab/bin"))


def usage() -> None:
    print(
        "Usage: python run_minibwa_hic_pipeline.py <hic_file1> <hic_file2> <genome.fasta>",
        file=sys.stderr,
    )


def q(path_or_text: str | Path) -> str:
    return shlex.quote(str(path_or_text))


def run(cmd: str, step_name: str) -> None:
    print(f"\n===== {step_name} =====", flush=True)
    print(cmd, flush=True)
    subprocess.run(["bash", "-lc", cmd], check=True)


def require(path: Path, label: str) -> None:
    if not path.exists():
        raise SystemExit(f"ERROR: missing {label}: {path}")


def main() -> int:
    if len(sys.argv) != 4:
        usage()
        return 1

    r1 = Path(sys.argv[1])
    r2 = Path(sys.argv[2])
    genome = Path(sys.argv[3])

    require(r1, "Hi-C read1")
    require(r2, "Hi-C read2")
    require(genome, "genome FASTA")
    require(MINIBWA, "minibwa")
    require(JUICER_TOOLS, "juicer_tools.jar")

    env_path = f"{HIC_AB_BIN}:{os.environ.get('PATH', '')}"
    os.environ["PATH"] = env_path

    prefix = genome.name
    for suffix in (".fasta", ".fa", ".fna"):
        if prefix.endswith(suffix):
            prefix = prefix[: -len(suffix)]
            break
    prefix = f"{prefix}.minibwa"

    chrom_sizes = Path(f"{prefix}.chrom.sizes")
    genome_fai = Path(f"{genome}.fai")

    Path("logs").mkdir(exist_ok=True)
    Path("tmp").mkdir(exist_ok=True)

    print("Input:")
    print(f"  R1: {r1}")
    print(f"  R2: {r2}")
    print(f"  genome: {genome}")
    print(f"  prefix: {prefix}")
    print(f"  threads: {THREADS}")
    print(f"  minibwa: {MINIBWA}")
    print(f"  juicer_tools: {JUICER_TOOLS}")

    run("command -v samtools && command -v pairtools && command -v pairix && command -v cooler && command -v java", "00 check tools")

    if not genome_fai.exists():
        run(f"samtools faidx {q(genome)}", "00 create FASTA index")
    run(f"cut -f1,2 {q(genome_fai)} > {q(chrom_sizes)} && ls -lh {q(chrom_sizes)}", "00 make chrom.sizes")

    if not Path(f"{genome}.mbw").exists() or not Path(f"{genome}.l2b").exists():
        run(f"{q(MINIBWA)} index {q(genome)} > logs/{q(prefix)}.minibwa.index.log 2>&1", "01 build minibwa index")
    else:
        print("\n===== 01 build minibwa index =====")
        print("minibwa index exists. Skip.")

    sorted_bam = Path(f"{prefix}.hic.sorted.bam")
    if not sorted_bam.exists():
        run(
            f"{q(MINIBWA)} mem --hic -t {THREADS} {q(genome)} {q(r1)} {q(r2)} "
            f"2> logs/{q(prefix)}.map.log "
            f"| samtools view -@ {THREADS} -bS - "
            f"| samtools sort -@ {THREADS} -o {q(sorted_bam)} - && "
            f"samtools index {q(sorted_bam)} && "
            f"samtools flagstat {q(sorted_bam)} > {q(prefix)}.hic.sorted.flagstat.txt",
            "02 minibwa mapping to coordinate-sorted BAM",
        )
    else:
        print("\n===== 02 minibwa mapping to coordinate-sorted BAM =====")
        print(f"{sorted_bam} exists. Skip.")

    namesort_bam = Path(f"{prefix}.hic.namesort.bam")
    if not namesort_bam.exists():
        run(
            f"samtools sort -n -@ {THREADS} -o {q(namesort_bam)} {q(sorted_bam)}",
            "03 name-sort BAM",
        )
    else:
        print("\n===== 03 name-sort BAM =====")
        print(f"{namesort_bam} exists. Skip.")

    dedup_pairsam = Path(f"{prefix}.dedup.pairsam.gz")
    if not dedup_pairsam.exists():
        run(
            f"samtools view -h {q(namesort_bam)} "
            f"| pairtools parse "
            f"--min-mapq 30 "
            f"--walks-policy 5unique "
            f"--max-inter-align-gap 30 "
            f"--chroms-path {q(chrom_sizes)} "
            f"--assembly {q(prefix)} "
            f"--output {q(prefix)}.pairsam.gz "
            f"--nproc-in {THREADS} "
            f"--nproc-out {THREADS} && "
            f"pairtools sort "
            f"--nproc {THREADS} "
            f"--tmpdir ./tmp "
            f"--output {q(prefix)}.sorted.pairsam.gz "
            f"{q(prefix)}.pairsam.gz && "
            f"pairtools dedup "
            f"--mark-dups "
            f"--output-stats {q(prefix)}.dedup.stats.txt "
            f"--output {q(dedup_pairsam)} "
            f"{q(prefix)}.sorted.pairsam.gz",
            "04 pairtools parse / sort / dedup",
        )
    else:
        print("\n===== 04 pairtools parse / sort / dedup =====")
        print(f"{dedup_pairsam} exists. Skip.")

    pairs_gz = Path(f"{prefix}.pairs.gz")
    if not pairs_gz.exists():
        run(
            f"pairtools split --output-pairs {q(pairs_gz)} --output-sam /dev/null {q(dedup_pairsam)}",
            "05 make pairs.gz",
        )
    else:
        print("\n===== 05 make pairs.gz =====")
        print(f"{pairs_gz} exists. Skip.")
    if not Path(f"{pairs_gz}.px2").exists():
        run(f"pairix {q(pairs_gz)}", "05 pairix index pairs.gz")

    cool = Path(f"{prefix}.10k.cool")
    if not cool.exists():
        run(
            f"cooler cload pairs "
            f"-c1 2 -p1 3 -c2 4 -p2 5 "
            f"{q(str(chrom_sizes) + ':10000')} "
            f"{q(pairs_gz)} "
            f"{q(cool)} && "
            f"cooler balance "
            f"--ignore-diags 2 "
            f"--mad-max 5 "
            f"--min-nnz 10 "
            f"-p {THREADS} "
            f"{q(cool)}",
            "06 make and balance 10k .cool",
        )
    else:
        print("\n===== 06 make and balance 10k .cool =====")
        print(f"{cool} exists. Skip.")
    run(f"cooler info {q(cool)} > {q(prefix)}.10k.cool.info.txt", "06 write cooler info")

    juicer_short = Path(f"{prefix}.juicer_short.txt")
    hic = Path(f"{prefix}.hic")
    if not juicer_short.exists():
        awk_cmd = (
            "awk 'BEGIN{OFS=\"\\t\"} /^#/ {next} "
            "{s1=($6==\"+\" ? 0 : 1); s2=($7==\"+\" ? 0 : 1); print s1,$2,$3,0,s2,$4,$5,1}'"
        )
        run(f"zcat {q(pairs_gz)} | {awk_cmd} > {q(juicer_short)}", "07 convert pairs.gz to Juicer short format")
    else:
        print("\n===== 07 convert pairs.gz to Juicer short format =====")
        print(f"{juicer_short} exists. Skip.")

    if not hic.exists():
        run(
            f"java -Xmx80g -jar {q(JUICER_TOOLS)} pre "
            f"-r 10000,20000,50000,100000,200000,500000,1000000 "
            f"{q(juicer_short)} "
            f"{q(hic)} "
            f"{q(chrom_sizes)} "
            f"> logs/{q(prefix)}.make_hic.log 2>&1",
            "07 make .hic for Juicebox",
        )
    else:
        print("\n===== 07 make .hic for Juicebox =====")
        print(f"{hic} exists. Skip.")

    print("\nDone.")
    print(f"Output cool: {cool}")
    print(f"Output hic:  {hic}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
