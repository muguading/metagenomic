#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(file2meco)
  library(microeco)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/run_microeco_beta.R <feature_table.qza> <taxonomy.qza>",
      "<sample_metadata.tsv> <group_column> <output_dir>",
      "[measure=bray]"
    )
  )
}

feature_table <- normalizePath(args[[1]], mustWork = TRUE)
taxonomy_table <- normalizePath(args[[2]], mustWork = TRUE)
sample_table <- normalizePath(args[[3]], mustWork = TRUE)
group_column <- args[[4]]
output_dir <- normalizePath(args[[5]], mustWork = FALSE)
measure <- if (length(args) >= 6) args[[6]] else "bray"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dataset <- qiime2meco(
  feature_table = feature_table,
  sample_table = sample_table,
  taxonomy_table = taxonomy_table
)

if (!group_column %in% colnames(dataset$sample_table)) {
  stop(sprintf("Group column not found in metadata: %s", group_column))
}

dataset$filter_pollution(taxa = c("mitochondria", "chloroplast"))
dataset$tidy_dataset(main_data = TRUE)
dataset$cal_betadiv(method = measure)

sample_counts <- as.data.frame(table(dataset$sample_table[[group_column]]), stringsAsFactors = FALSE)
colnames(sample_counts) <- c(group_column, "sample_count")
sample_counts <- sample_counts[sample_counts$sample_count > 0, , drop = FALSE]

beta_obj <- trans_beta$new(
  dataset = dataset,
  measure = measure,
  group = group_column
)

nmds_obj <- trans_beta$new(
  dataset = dataset,
  measure = measure,
  group = group_column
)

beta_obj$cal_ordination(method = "PCoA")
nmds_obj$cal_ordination(method = "NMDS")
beta_obj$cal_manova(manova_all = TRUE)
beta_obj$cal_anosim()
beta_obj$cal_betadisper()
beta_obj$cal_group_distance(within_group = TRUE)
beta_obj$cal_group_distance_diff()

group_levels <- unique(as.character(dataset$sample_table[[group_column]]))
palette_values <- scales::hue_pal()(length(group_levels))
names(palette_values) <- group_levels

pcoa_plot <- beta_obj$plot_ordination(
  plot_color = group_column,
  plot_type = "point",
  color_values = palette_values,
  point_size = 2.4,
  point_alpha = 0.85
)

nmds_plot <- nmds_obj$plot_ordination(
  plot_color = group_column,
  plot_type = "point",
  color_values = palette_values,
  point_size = 2.4,
  point_alpha = 0.85,
  NMDS_stress_text_prefix = "Stress=="
)

distance_plot <- beta_obj$plot_group_distance(plot_group_order = sort(group_levels))

write.table(
  dataset$sample_table,
  file = file.path(output_dir, "sample_table_used.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  sample_counts,
  file = file.path(output_dir, "group_sample_counts.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  as.matrix(dataset$beta_diversity[[measure]]),
  file = file.path(output_dir, sprintf("%s_distance_matrix.tsv", measure)),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  beta_obj$res_ordination$scores,
  file = file.path(output_dir, "pcoa_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  nmds_obj$res_ordination$scores,
  file = file.path(output_dir, "nmds_scores.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  beta_obj$res_manova,
  file = file.path(output_dir, "permanova.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  beta_obj$res_anosim,
  file = file.path(output_dir, "anosim.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  beta_obj$res_group_distance,
  file = file.path(output_dir, "within_group_distances.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)
write.table(
  beta_obj$res_group_distance_diff,
  file = file.path(output_dir, "within_group_distance_stats.tsv"),
  sep = "\t",
  quote = FALSE,
  col.names = NA
)

capture.output(beta_obj$res_betadisper, file = file.path(output_dir, "betadisper.txt"))
saveRDS(dataset, file = file.path(output_dir, "microeco_dataset.rds"))
saveRDS(beta_obj, file = file.path(output_dir, "microeco_trans_beta.rds"))
saveRDS(nmds_obj, file = file.path(output_dir, "microeco_trans_beta_nmds.rds"))

ggsave(
  filename = file.path(output_dir, "pcoa_plot.png"),
  plot = pcoa_plot,
  width = 9,
  height = 7,
  dpi = 300
)
ggsave(
  filename = file.path(output_dir, "pcoa_plot.pdf"),
  plot = pcoa_plot,
  width = 9,
  height = 7
)
ggsave(
  filename = file.path(output_dir, "nmds_plot.png"),
  plot = nmds_plot,
  width = 9,
  height = 7,
  dpi = 300
)
ggsave(
  filename = file.path(output_dir, "nmds_plot.pdf"),
  plot = nmds_plot,
  width = 9,
  height = 7
)
ggsave(
  filename = file.path(output_dir, "within_group_distance_plot.png"),
  plot = distance_plot,
  width = 10,
  height = 7,
  dpi = 300
)
ggsave(
  filename = file.path(output_dir, "within_group_distance_plot.pdf"),
  plot = distance_plot,
  width = 10,
  height = 7
)

summary_lines <- c(
  sprintf("feature_table: %s", feature_table),
  sprintf("taxonomy_table: %s", taxonomy_table),
  sprintf("sample_table: %s", sample_table),
  sprintf("group_column: %s", group_column),
  sprintf("measure: %s", measure),
  sprintf("samples_used: %s", nrow(dataset$sample_table)),
  sprintf("features_used: %s", nrow(dataset$otu_table)),
  sprintf("taxa_rows_used: %s", nrow(dataset$tax_table)),
  sprintf("group_levels: %s", paste(sort(group_levels), collapse = ", ")),
  sprintf("nmds_stress: %s", as.character(nmds_obj$res_ordination$model$stress))
)

writeLines(summary_lines, con = file.path(output_dir, "run_summary.txt"))
