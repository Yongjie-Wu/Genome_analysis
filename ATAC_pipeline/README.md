# Single-haplotype ATAC-seq pipeline

This pipeline analyzes one ATAC-seq dataset against one reference genome, for example only `hap1` or only `hap2`. If you want to analyze another haplotype, run the same pipeline again with another genome FASTA/GFF and another project name.

## Step scripts

Run any script without arguments to see its usage.

```bash
bash 00.clean.sh <sample_name> <R1.fq.gz> <R2.fq.gz> [project_dir]
bash 01.index.sh <genome.fa> [project_dir]
bash 02.map.sh <sample_name> <clean_R1.fq.gz> <clean_R2.fq.gz> <genome.fa> [project_dir]
bash 03.filter_markdup_merge.sh <project_name> <project_dir> <sample1> <sample2> [sample3 ...]
bash 04.macs2_callpeak.sh <project_name> <genome_size> [project_dir]
bash 05.make_bigwig.sh <project_name> [project_dir]
bash 06.qc_peak_signal.sh <project_name> <project_dir> <sample1> <sample2> [sample3 ...]
bash 07.peak_gene_annotation.sh <project_name> <gene.gff3> [project_dir]
bash 08.atac_window_features.sh <project_name> <genome.fa> <gene.gff3> <centromere_range.txt> <SatDNA.bed> <RPKM.tsv> [project_dir]
bash 09.summary.sh <project_name> [project_dir]
```

## Run all

`samples.tsv` format:

```text
sample	R1	R2
CRRxxxx	CRRxxxx_r1.fq.gz	CRRxxxx_r2.fq.gz
```

Run standard ATAC analysis:

```bash
bash run_all.sh duli_hap1 samples.tsv hap1/Pyrbe_duli_hap1.fa hap1/Pyrbe_duli_hap1.gff3 5e8 .
```

Run with ATAC/SatDNA/expression window analysis:

```bash
bash run_all.sh duli_hap1 samples.tsv hap1/Pyrbe_duli_hap1.fa hap1/Pyrbe_duli_hap1.gff3 5e8 . \
  --step08-extra Pyrbe_duli.hap1_centromere_range.txt Pyrbe_duli.hap1.SatDNA.bed Pyrbe_duli_hap1.RPKM
```

If clean reads already exist:

```bash
bash run_all.sh duli_hap1 samples.clean.tsv hap1/Pyrbe_duli_hap1.fa hap1/Pyrbe_duli_hap1.gff3 5e8 . --skip-clean
```

## `genome_size`

`genome_size` is the effective genome size used by MACS2, passed to `macs2 callpeak -g`. For pear haplotype genomes, `5e8` is a reasonable default.

## Conda environments

Default environment variables:

```bash
WORK_ENV=work
MACS2_ENV=macs2_clean
DEEPTOOLS_ENV=deeptools_local
PY_ENV=base
```
