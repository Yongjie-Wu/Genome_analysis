# Hi-C Triangle Plus Multi-Track Plot

This workflow makes a figure with a Hi-C triangular contact map on top and one-chromosome genomic feature tracks below.

## Step 1. Prepare Tracks

```bash
bash prepare_hic_triangle_tracks.sh <sample_name> <chromosome> <genome.fa.fai> <centromere_range.txt> <Gypsy.bed> <Copia.bed> <gene.gff3|Gene.bed> <SatDNA.bed> <CpG.cov.bed> <CHG.cov.bed> <CHH.cov.bed> <output_dir> [options]
```

Optional:

```bash
--rpkm <RPKM.tsv>
--gypsy-ltr <Gypsy_LTR.bed>
--gypsy-int <Gypsy_INT.bed>
--window 30000
```

## Step 2. Plot

```bash
bash plot_hic_triangle_tracks.sh <cool_file> <chromosome> <centromere_range.txt> <Copia_track.txt> <Gypsy_track.txt> <SatDNA_track.txt> <Gene_track.txt> <LeafFruit_expression.tsv|none> <CG_track.txt> <CHG_track.txt> <CHH_track.txt> <output_prefix> [options]
```

Optional:

```bash
--gypsy-ltr <Gypsy_LTR.bed>
--gypsy-int <Gypsy_INT.bed>
--dpi 600
--no-pdf
--hic-bin-step 2
--font-file /path/to/Arial.ttf
```

## Example

```bash
bash prepare_hic_triangle_tracks.sh Pyrbe_duli_hap1 Chr02A \
  /path/to/Pyrbe_duli_hap1.fa.fai \
  /path/to/Pyrbe_duli.hap1_centromere_range.txt \
  /path/to/Gypsy.bed \
  /path/to/Copia.bed \
  /path/to/Pyrbe_duli_hap1.gff3 \
  /path/to/Pyrbe_duli.hap1.SatDNA.bed \
  /path/to/hap1.CpG.cov.bed \
  /path/to/hap1.CHG.cov.bed \
  /path/to/hap1.CHH.cov.bed \
  hic_tracks_input \
  --rpkm /path/to/Pyrbe_duli_hap1.RPKM \
  --gypsy-ltr /path/to/Chr02A_Gypsy_LTR.bed \
  --gypsy-int /path/to/Chr02A_Gypsy_INT.bed

bash plot_hic_triangle_tracks.sh \
  /path/to/Pyrbe_duli_hap1.10k.cool \
  Chr02A \
  /path/to/Pyrbe_duli.hap1_centromere_range.txt \
  hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Copia.30000bp.txt \
  hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Gypsy.30000bp.txt \
  hic_tracks_input/Pyrbe_duli_hap1.Chr02A.SatDNA.30000bp.txt \
  hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Gene.30000bp.txt \
  hic_tracks_input/Pyrbe_duli_hap1.Chr02A.leaf_fruit_expression.30000bp.tsv \
  hic_tracks_input/Pyrbe_duli_hap1.Chr02A.CG.30000bp.txt \
  hic_tracks_input/Pyrbe_duli_hap1.Chr02A.CHG.30000bp.txt \
  hic_tracks_input/Pyrbe_duli_hap1.Chr02A.CHH.30000bp.txt \
  output/Chr02A_hic_triangle_plus_multi_tracks \
  --gypsy-ltr hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Gypsy_LTR.bed \
  --gypsy-int hic_tracks_input/Pyrbe_duli_hap1.Chr02A.Gypsy_INT.bed
```
