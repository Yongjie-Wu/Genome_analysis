#!/usr/bin/env Rscript

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript plot_chr_tracks.R <sample_name> <chr_id> <window_size> <project_dir>\n")
  quit(status = 1)
}

sample_name <- args[1]
chr_id <- args[2]
win <- args[3]
project_dir <- args[4]

read_track <- function(path, value_name) {
  if (!file.exists(path)) stop("missing file: ", path)
  x <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(x)[1:4] <- c("chr", "start", "end", value_name)
  x <- x[x$chr == chr_id, c("chr", "start", "end", value_name)]
  x$mid <- (x$start + x$end) / 2 / 1e6
  x
}

read_features <- function(path) {
  if (!file.exists(path)) stop("missing file: ", path)
  x <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  x <- x[x$chr == chr_id, ]
  x$mid <- (x$start + x$end) / 2 / 1e6
  x
}

cg <- read_track(file.path(project_dir, "03.windows", paste0(sample_name, ".CG.", win, ".bed")), "value")
chg <- read_track(file.path(project_dir, "03.windows", paste0(sample_name, ".CHG.", win, ".bed")), "value")
chh <- read_track(file.path(project_dir, "03.windows", paste0(sample_name, ".CHH.", win, ".bed")), "value")
features <- read_features(file.path(project_dir, "04.features", paste0(sample_name, ".feature_density.", win, ".tsv")))

xmax <- max(c(cg$end, chg$end, chh$end, features$end), na.rm = TRUE) / 1e6
cen <- features[features$in_centromere == 1, ]
cen_start <- if (nrow(cen) > 0) min(cen$start) / 1e6 else NA
cen_end <- if (nrow(cen) > 0) max(cen$end) / 1e6 else NA

shade <- function() {
  if (is.na(cen_start)) {
    NULL
  } else {
    annotate("rect", xmin = cen_start, xmax = cen_end, ymin = -Inf, ymax = Inf,
             fill = "grey85", alpha = 0.5)
  }
}

theme_track <- theme_classic(base_size = 9) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.ticks.y = element_line(linewidth = 0.2),
    axis.line = element_line(linewidth = 0.25),
    plot.margin = margin(1, 4, 1, 4)
  )

bar_track <- function(df, y, label, color, ymax = NULL) {
  p <- ggplot(df, aes(x = mid, y = .data[[y]])) +
    shade() +
    geom_col(width = as.numeric(win) / 1e6, fill = color, linewidth = 0) +
    scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
    labs(x = NULL, y = label) +
    theme_track
  if (!is.null(ymax)) p <- p + coord_cartesian(ylim = c(0, ymax))
  p
}

p_gene <- bar_track(features, "gene_density", "Gene", "#228B22", 1)
p_sat <- bar_track(features, "satdna_density", "SatDNA", "#72329D", 1)
p_copia <- bar_track(features, "copia_density", "Copia", "#FF1418", 1)
p_gypsy <- bar_track(features, "gypsy_density", "Gypsy", "#00AEEE", 1)
p_cg <- bar_track(cg, "value", "CG", "#C80F12", 1)
p_chg <- bar_track(chg, "value", "CHG", "#0012FF", 1)
p_chh <- bar_track(chh, "value", "CHH", "#F59646", 0.3)

plot <- p_copia / p_gypsy / p_sat / p_gene / p_cg / p_chg / p_chh +
  plot_layout(heights = rep(1, 7)) &
  theme(axis.title.x = element_blank())

plot <- plot & labs(x = paste0("Position on ", chr_id, " (Mb)"))

outdir <- file.path(project_dir, "05.plots")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
pdf_file <- file.path(outdir, paste0(sample_name, ".", chr_id, ".methylation_tracks.pdf"))
png_file <- file.path(outdir, paste0(sample_name, ".", chr_id, ".methylation_tracks.png"))

ggsave(pdf_file, plot, width = 7, height = 8)
ggsave(png_file, plot, width = 7, height = 8, dpi = 300)
message("Wrote ", pdf_file)
message("Wrote ", png_file)
