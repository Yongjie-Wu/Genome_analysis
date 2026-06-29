# Genome feature circos plot

This script draws a 17-chromosome multi-track circos plot from prepared window tracks.

## Required files

1. Genome FASTA index: `genome.fa.fai`
2. Centromere range file:

```text
Chr01A  start  end
Chr02A  start  end
```

3. Genome-wide track files. Each track file has four columns:

```text
chr  start  end  value
```

Required tracks:

- Gypsy density
- Copia density
- Gene density
- SatDNA density
- CG methylation
- CHG methylation
- CHH methylation

Expression and ATAC tracks are optional. Expression file must have:

```text
chr  start  end  leaf  fruit
```

Tracks drawn:

- outer chromosome ring
- centromere positions, drawn as black blocks on the chromosome ring
- Copia density
- Gypsy density
- SatDNA density
- ATAC signal, if `--atac` is provided
- Gene density
- leaf expression, if expression files exist
- fruit expression, if expression files exist
- CG methylation
- CHG methylation
- CHH methylation

## Run

First prepare genome-wide window tables from basic files:

```bash
bash prepare_circos_inputs.sh <sample_name> <genome.fa.fai> <centromere_range.txt> <Gypsy.bed> <Copia.bed> <gene.gff3|Gene.bed> <SatDNA.bed> <CpG.cov.bed> <CHG.cov.bed> <CHH.cov.bed> <output_dir> [options]
```

Then draw the circos plot:

```bash
bash run_genome_feature_circos.sh <genome.fa.fai> <centromere_range.txt> <Gypsy_30kb.txt> <Copia_30kb.txt> <Gene_30kb.txt> <SatDNA_30kb.txt> <CG_30kb.txt> <CHG_30kb.txt> <CHH_30kb.txt> <output_prefix> [options]
```

Example:

```bash
bash prepare_circos_inputs.sh Pyrbe_duli_hap1 \
  Pyrbe_duli_hap1.fa.fai \
  Pyrbe_duli.hap1_centromere_range.txt \
  Gypsy.bed \
  Copia.bed \
  Pyrbe_duli_hap1.gff3 \
  Pyrbe_duli.hap1.SatDNA.bed \
  hap1.CpG.cov.bed \
  hap1.CHG.cov.bed \
  hap1.CHH.cov.bed \
  circos_input \
  --rpkm Pyrbe_duli_hap1.RPKM \
  --atac-bam /path/to/duli.hap1.merge.rmdup.bam

bash run_genome_feature_circos.sh \
  Pyrbe_duli_hap1.fa.fai \
  Pyrbe_duli.hap1_centromere_range.txt \
  circos_input/Pyrbe_duli_hap1.Gypsy.30000bp.txt \
  circos_input/Pyrbe_duli_hap1.Copia.30000bp.txt \
  circos_input/Pyrbe_duli_hap1.Gene.30000bp.txt \
  circos_input/Pyrbe_duli_hap1.SatDNA.30000bp.txt \
  circos_input/Pyrbe_duli_hap1.CG.30000bp.txt \
  circos_input/Pyrbe_duli_hap1.CHG.30000bp.txt \
  circos_input/Pyrbe_duli_hap1.CHH.30000bp.txt \
  output/Pyrbe_duli_hap1_17chr_circos \
  --atac circos_input/Pyrbe_duli_hap1.ATAC.30000bp.tsv \
  --expression circos_input/Pyrbe_duli_hap1.leaf_fruit_expression.30000bp.tsv
```

Example:

```bash
bash run_genome_feature_circos.sh \
  /path/to/project_dir \
  /path/to/genome.fa.fai \
  /path/to/centromere_range.txt \
  /path/to/output/sample_17chr_circos \
  --atac /path/to/ATAC_window.tsv \
  --title "sample hap1 17-chromosome multi-track circos"
```

Outputs:

```text
sample_17chr_circos.png
sample_17chr_circos.pdf
sample_17chr_circos.svg
```
