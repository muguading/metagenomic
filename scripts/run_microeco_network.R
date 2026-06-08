#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(file2meco)
  library(microeco)
  library(igraph)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/run_microeco_network.R <feature_table.qza> <taxonomy.qza>",
      "<sample_metadata.tsv> <group_column> <output_dir>",
      "[taxa_level=Genus] [cor_method=spearman] [cor_cut=0.6] [p_thres=0.05] [filter_thres=0.0005]"
    )
  )
}

feature_table <- normalizePath(args[[1]], mustWork = TRUE)
taxonomy_table <- normalizePath(args[[2]], mustWork = TRUE)
sample_table <- normalizePath(args[[3]], mustWork = TRUE)
group_column <- args[[4]]
output_dir <- normalizePath(args[[5]], mustWork = FALSE)
taxa_level <- if (length(args) >= 6) args[[6]] else "Genus"
cor_method <- if (length(args) >= 7) args[[7]] else "spearman"
cor_cut <- if (length(args) >= 8) as.numeric(args[[8]]) else 0.6
p_thres <- if (length(args) >= 9) as.numeric(args[[9]]) else 0.05
filter_thres <- if (length(args) >= 10) as.numeric(args[[10]]) else 0.0005

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

normalize_taxa_level <- function(x) {
  text <- trimws(as.character(x))
  if (identical(text, "")) return("Genus")
  mapped <- switch(
    tolower(text),
    otu = "OTU",
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

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

write_empty_outputs <- function(reason = "network_not_available") {
  write.table(
    data.frame(metric = c("status", "reason"), value = c("empty", reason), stringsAsFactors = FALSE),
    file = file.path(output_dir, "network_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  write.table(data.frame(), file = file.path(output_dir, "node_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(), file = file.path(output_dir, "edge_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(), file = file.path(output_dir, "module_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(), file = file.path(output_dir, "role_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(data.frame(), file = file.path(output_dir, "eigen_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  writeLines(c(sprintf("status: empty"), sprintf("reason: %s", reason)), con = file.path(output_dir, "run_summary.txt"))
}

target_level <- normalize_taxa_level(taxa_level)

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

network <- trans_network$new(
  dataset = dataset,
  cor_method = cor_method,
  taxa_level = target_level,
  filter_thres = filter_thres
)

network$cal_network(
  COR_p_thres = p_thres,
  COR_cut = cor_cut,
  add_taxa_name = c("Phylum", "Class", "Order", "Family", "Genus"),
  delete_unlinked_nodes = TRUE
)

if (is.null(network$res_network) || igraph::gorder(network$res_network) < 2 || igraph::gsize(network$res_network) < 1) {
  write_empty_outputs("insufficient_edges_after_filtering")
  quit(save = "no", status = 0)
}

network$cal_module()
network$get_node_table(node_roles = TRUE)
network$get_edge_table()
network$cal_eigen()
network$cal_sum_links(taxa_level = "Phylum")

node_table <- network$res_node_table
edge_table <- network$res_edge_table

if (is.null(node_table) || nrow(node_table) < 1 || is.null(edge_table) || nrow(edge_table) < 1) {
  write_empty_outputs("empty_node_or_edge_table")
  quit(save = "no", status = 0)
}

node_table$name <- as.character(node_table$name)
edge_table$node1 <- as.character(edge_table$node1)
edge_table$node2 <- as.character(edge_table$node2)
edge_table$label <- as.character(edge_table$label)
edge_table$weight <- as.numeric(edge_table$weight)

graph <- network$res_network
membership <- if ("module" %in% colnames(node_table)) {
  node_table$module[match(igraph::V(graph)$name, node_table$name)]
} else {
  NULL
}

edge_weights <- edge_table$weight
positive_edges <- sum(edge_table$label == "+", na.rm = TRUE)
negative_edges <- sum(edge_table$label == "-", na.rm = TRUE)

module_summary <- do.call(
  rbind,
  lapply(split(node_table, node_table$module), function(df) {
    data.frame(
      module = as.character(df$module[[1]]),
      node_count = nrow(df),
      mean_degree = round(mean(as.numeric(df$degree), na.rm = TRUE), 4),
      max_degree = max(as.numeric(df$degree), na.rm = TRUE),
      total_abundance = round(sum(as.numeric(df$Abundance), na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  })
)
if (!is.null(module_summary) && nrow(module_summary) > 0) {
  module_summary <- module_summary[order(module_summary$node_count, decreasing = TRUE), , drop = FALSE]
} else {
  module_summary <- data.frame()
}

role_summary <- do.call(
  rbind,
  lapply(split(node_table, node_table$taxa_roles), function(df) {
    data.frame(
      taxa_role = as.character(df$taxa_roles[[1]]),
      node_count = nrow(df),
      mean_degree = round(mean(as.numeric(df$degree), na.rm = TRUE), 4),
      stringsAsFactors = FALSE
    )
  })
)
if (is.null(role_summary) || nrow(role_summary) < 1) {
  role_summary <- data.frame()
}

eigen_summary <- data.frame()
if (!is.null(network$res_eigen_expla)) {
  eigen_values <- network$res_eigen_expla
  if (is.null(dim(eigen_values))) {
    numeric_values <- suppressWarnings(as.numeric(gsub("%", "", as.character(eigen_values), fixed = TRUE)))
    eigen_summary <- data.frame(
      module = names(eigen_values) %||% seq_along(eigen_values),
      variance_explained = numeric_values,
      stringsAsFactors = FALSE
    )
  } else {
    eigen_summary <- as.data.frame(eigen_values, stringsAsFactors = FALSE)
    if (!"module" %in% colnames(eigen_summary)) {
      eigen_summary$module <- rownames(eigen_summary)
    }
  }
}

summary_table <- data.frame(
  metric = c(
    "status",
    "taxa_level",
    "cor_method",
    "cor_cut",
    "p_thres",
    "filter_thres",
    "node_count",
    "edge_count",
    "positive_edge_count",
    "negative_edge_count",
    "module_count",
    "avg_degree",
    "max_degree",
    "mean_abs_weight",
    "density",
    "transitivity",
    "modularity"
  ),
  value = c(
    "ready",
    target_level,
    cor_method,
    cor_cut,
    p_thres,
    filter_thres,
    nrow(node_table),
    nrow(edge_table),
    positive_edges,
    negative_edges,
    length(unique(node_table$module)),
    round(mean(as.numeric(node_table$degree), na.rm = TRUE), 4),
    max(as.numeric(node_table$degree), na.rm = TRUE),
    round(mean(abs(edge_weights), na.rm = TRUE), 4),
    round(igraph::edge_density(graph, loops = FALSE), 6),
    round(igraph::transitivity(graph, type = "globalundirected"), 6),
    if (!is.null(membership)) round(igraph::modularity(graph, membership = as.factor(membership)), 6) else NA_real_
  ),
  stringsAsFactors = FALSE
)

write.table(summary_table, file = file.path(output_dir, "network_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(node_table, file = file.path(output_dir, "node_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(edge_table, file = file.path(output_dir, "edge_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(module_summary, file = file.path(output_dir, "module_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(role_summary, file = file.path(output_dir, "role_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(eigen_summary, file = file.path(output_dir, "eigen_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

if (!is.null(network$res_sum_links_pos)) {
  write.table(as.data.frame(network$res_sum_links_pos), file = file.path(output_dir, "phylum_links_positive.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}
if (!is.null(network$res_sum_links_neg)) {
  write.table(as.data.frame(network$res_sum_links_neg), file = file.path(output_dir, "phylum_links_negative.tsv"), sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
}

summary_lines <- c(
  sprintf("status: ready"),
  sprintf("feature_table: %s", feature_table),
  sprintf("taxonomy_table: %s", taxonomy_table),
  sprintf("sample_table: %s", sample_table),
  sprintf("group_column: %s", group_column),
  sprintf("taxa_level: %s", target_level),
  sprintf("cor_method: %s", cor_method),
  sprintf("cor_cut: %s", cor_cut),
  sprintf("p_thres: %s", p_thres),
  sprintf("filter_thres: %s", filter_thres),
  sprintf("node_count: %s", nrow(node_table)),
  sprintf("edge_count: %s", nrow(edge_table)),
  sprintf("module_count: %s", length(unique(node_table$module))),
  sprintf("positive_edge_count: %s", positive_edges),
  sprintf("negative_edge_count: %s", negative_edges)
)

writeLines(summary_lines, con = file.path(output_dir, "run_summary.txt"))
