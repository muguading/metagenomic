#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(file2meco)
  library(microeco)
  library(ggplot2)
  library(randomForest)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/run_microeco_biomarker.R <feature_table.qza> <taxonomy.qza>",
      "<sample_metadata.tsv> <group_column> <output_dir>",
      "[taxa_level=Genus] [rf_model=random_forest]"
    )
  )
}

feature_table <- normalizePath(args[[1]], mustWork = TRUE)
taxonomy_table <- normalizePath(args[[2]], mustWork = TRUE)
sample_table <- normalizePath(args[[3]], mustWork = TRUE)
group_column <- args[[4]]
output_dir <- normalizePath(args[[5]], mustWork = FALSE)
taxa_level <- if (length(args) >= 6) args[[6]] else "Genus"
rf_model <- if (length(args) >= 7) args[[7]] else "random_forest"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalize_taxa_level <- function(x) {
  text <- trimws(as.character(x))
  if (identical(text, "")) return("Genus")
  mapped <- switch(
    tolower(text),
    phylum = "Phylum",
    class = "Class",
    order = "Order",
    family = "Family",
    genus = "Genus",
    species = "Species",
    text
  )
  mapped
}

clean_taxon_label <- function(x) {
  text <- trimws(as.character(x))
  text <- sub("^.*\\|", "", text)
  text <- sub("^[a-z]__", "", text)
  text <- gsub("_", " ", text, fixed = TRUE)
  text <- trimws(text)
  if (identical(text, "")) {
    return("未注释")
  }
  text
}

is_annotated_taxon <- function(x) {
  text <- tolower(trimws(as.character(x)))
  !(
    identical(text, "") ||
      text %in% c("未注释", "未分类", "其他", "unclassified", "unassigned", "unknown", "norank", "uncultured", "ambiguous taxa") ||
      startsWith(text, "unclassified ") ||
      startsWith(text, "unknown ")
  )
}

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
dataset$cal_abund()

target_level <- normalize_taxa_level(taxa_level)

lefse_obj <- trans_diff$new(
  dataset = dataset,
  method = "lefse",
  group = group_column,
  taxa_level = target_level,
  remove_unknown = TRUE
)

lefse_table <- lefse_obj$res_diff
write.table(
  lefse_table,
  file = file.path(output_dir, "lefse_diff.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (nrow(lefse_table) > 0) {
  top_n <- min(20, nrow(lefse_table))
  lefse_bar <- lefse_obj$plot_diff_bar(use_number = 1:top_n)
  ggsave(
    filename = file.path(output_dir, "lefse_barplot.png"),
    plot = lefse_bar,
    width = 10,
    height = 7,
    dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, "lefse_barplot.pdf"),
    plot = lefse_bar,
    width = 10,
    height = 7
  )
}

abundance_table <- dataset$taxa_abund[[target_level]]
if (is.null(abundance_table)) {
  stop(sprintf("microeco taxa_abund does not contain level: %s", target_level))
}

sample_ids <- intersect(colnames(abundance_table), rownames(dataset$sample_table))
rf_groups <- as.character(dataset$sample_table[sample_ids, group_column, drop = TRUE])
rf_keep <- !is.na(rf_groups) & trimws(rf_groups) != ""
sample_ids <- sample_ids[rf_keep]
rf_groups <- rf_groups[rf_keep]

group_counts <- table(rf_groups)
valid_groups <- names(group_counts[group_counts >= 2])
group_keep <- rf_groups %in% valid_groups
sample_ids <- sample_ids[group_keep]
rf_groups <- rf_groups[group_keep]

rf_matrix <- as.matrix(abundance_table[, sample_ids, drop = FALSE])
taxa_labels <- rownames(rf_matrix)
clean_labels <- vapply(taxa_labels, clean_taxon_label, character(1))
annotated_keep <- vapply(clean_labels, is_annotated_taxon, logical(1))
rf_matrix <- rf_matrix[annotated_keep, , drop = FALSE]
clean_labels <- clean_labels[annotated_keep]

if (nrow(rf_matrix) > 1 && length(unique(rf_groups)) >= 2) {
  collapsed <- rowsum(rf_matrix, group = clean_labels, reorder = FALSE)
  rf_input <- as.data.frame(t(collapsed))
  rf_fit <- randomForest(
    x = rf_input,
    y = factor(rf_groups),
    importance = TRUE,
    ntree = 500
  )
  importance_matrix <- randomForest::importance(rf_fit)
  if (is.matrix(importance_matrix)) {
    if ("MeanDecreaseAccuracy" %in% colnames(importance_matrix)) {
      importance_values <- importance_matrix[, "MeanDecreaseAccuracy"]
    } else if ("MeanDecreaseGini" %in% colnames(importance_matrix)) {
      importance_values <- importance_matrix[, "MeanDecreaseGini"]
    } else {
      importance_values <- importance_matrix[, ncol(importance_matrix)]
    }
  } else {
    importance_values <- importance_matrix
  }
  rf_table <- data.frame(
    Taxa = rownames(importance_matrix),
    Importance = as.numeric(importance_values),
    Method = rf_model,
    stringsAsFactors = FALSE
  )
  rf_table <- rf_table[order(rf_table$Importance, decreasing = TRUE), , drop = FALSE]
  write.table(
    rf_table,
    file = file.path(output_dir, "rf_importance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  if (nrow(rf_table) > 0) {
    top_rf <- head(rf_table, 20)
    rf_plot <- ggplot(top_rf, aes(x = reorder(Taxa, Importance), y = Importance)) +
      geom_col(fill = "#6f8f84", width = 0.72) +
      coord_flip() +
      theme_minimal(base_size = 12) +
      labs(x = NULL, y = "Feature importance", title = sprintf("Random Forest Biomarker (%s)", target_level))
    ggsave(
      filename = file.path(output_dir, "rf_importance.png"),
      plot = rf_plot,
      width = 10,
      height = 7,
      dpi = 300
    )
    ggsave(
      filename = file.path(output_dir, "rf_importance.pdf"),
      plot = rf_plot,
      width = 10,
      height = 7
    )
  }
} else {
  write.table(
    data.frame(Taxa = character(), Importance = numeric(), Method = character(), stringsAsFactors = FALSE),
    file = file.path(output_dir, "rf_importance.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

summary_lines <- c(
  sprintf("feature_table: %s", feature_table),
  sprintf("taxonomy_table: %s", taxonomy_table),
  sprintf("sample_table: %s", sample_table),
  sprintf("group_column: %s", group_column),
  sprintf("taxa_level: %s", target_level),
  sprintf("rf_model: %s", rf_model),
  sprintf("samples_used: %s", length(sample_ids)),
  sprintf("groups_used: %s", paste(sort(unique(rf_groups)), collapse = ", ")),
  sprintf("lefse_hits: %s", nrow(lefse_table))
)

writeLines(summary_lines, con = file.path(output_dir, "run_summary.txt"))
