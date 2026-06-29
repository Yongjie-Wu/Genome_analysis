# Single-haplotype methylation pipeline

This pipeline processes one haplotype/reference at a time.

Required input files:

- PacBio HiFi modification BAM, for example `mod.bam`
- one haplotype genome FASTA, for example `hap1/*.fa`
- one haplotype gene annotation GFF3
- centromere range file: `chr start end`
- SatDNA BED
- optional Copia/Gypsy BED files, or use `none`

Main command:

```bash
bash run_all.sh <sample_name> <mod.bam> <genome.fa> <gene.gff3> <centromere_range.txt> <SatDNA.bed> <Copia.bed|none> <Gypsy.bed|none> <genome_ref_prefix> [project_dir] [options]
```

Example:

```bash
bash run_all.sh dianli_hap1 \
  /storage3/wuyongjie/Project/centromere/jiajihua/dianli/mod.bam \
  /storage3/wuyongjie/Project/centromere/jiajihua/dianli/hap1/Pyrps_dianli_hap1.fa \
  /storage3/wuyongjie/Project/centromere/jiajihua/dianli/hap1/Pyrps_dianli_hap1.gff3 \
  /storage3/wuyongjie/Project/centromere/jiajihua/dianli/hap1/0.Pyrps_dianli.hap1.FINAL.fixTelo_range.txt \
  /storage3/wuyongjie/Project/centromere/jiajihua/dianli/hap1/Pyrps_dianli.hap1.SatDNA.bed \
  /storage3/wuyongjie/Project/centromere/jiajihua/dianli/hap1/Copia.bed \
  /storage3/wuyongjie/Project/centromere/jiajihua/dianli/hap1/Gypsy.bed \
  dianli_hap1 \
  /storage3/wuyongjie/Project/centromere/jiajihua/dianli_methylation_hap1 \
  --threads 48 --window 30000 --plot-chr Chr02A
```

Step scripts:

```bash
bash 00.index_reference.sh <genome.fa> <ref_prefix> [project_dir]
bash 01.align_modbam.sh <sample_name> <mod.bam> <genome.fa> <ref.mmi> [project_dir] [threads]
bash 02.hifimeth_pileup.sh <sample_name> <genome.fa> <aligned.pbmm2.bam> [project_dir]
bash 03.window_methylation.sh <sample_name> <genome.fa> <CpG.cov.bed> <CHG.cov.bed> <CHH.cov.bed> [window_size] [project_dir]
bash 04.feature_density.sh <sample_name> <genome.fa> <gene.gff3> <centromere_range.txt> <SatDNA.bed> <Copia.bed|none> <Gypsy.bed|none> [window_size] [project_dir]
bash 05.plot_chr_tracks.sh <sample_name> <chr_id> [window_size] [project_dir]
bash 06.summary.sh <sample_name> [project_dir]
```

Outputs:

- `00.index`: FASTA / pbmm2 index files
- `01.align`: sorted pbmm2 BAM
- `02.methylation`: hifimeth CpG/CHG/CHH cov.bed
- `03.windows`: weighted methylation level per window
- `04.features`: Gene/SatDNA/Copia/Gypsy density and centromere flag per window
- `05.plots`: single-chromosome track plots
- `06.summary`: output file summary
