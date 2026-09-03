#' ##########################################################################
#' guide_representation.R
#'
#' Generates QC and representation plots (ECDF, Lorenz/Gini, Density, Dropout)
#' from MAGeCK count tables.
#'
#' USAGE:
#'    Rscript Scripts/utils/guide_representation.R [path_to_count_table.txt]
#'
#' If no path is provided, automatically processes all *.count.txt files
#' found in Analysis_Data/counts/.
#' ##########################################################################

options(stringsAsFactors = FALSE)

# Load core R packages
suppressMessages({
  library(dplyr, quietly = TRUE)
  library(tidyr, quietly = TRUE)
  library(ggplot2, quietly = TRUE)
})

# Custom Gini index calculation (avoid third-party CRAN dependency)
calc_gini <- function(x) {
  x <- sort(x[x > 0])
  n <- length(x)
  if (n == 0 || sum(x) == 0) {
    return(NA)
  }
  sum((2 * seq_len(n) - n - 1) * x) / (n * sum(x))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  count_files <- args
} else {
  count_files <- list.files("Analysis_Data/counts", pattern = "\\.count\\.txt$", full.names = TRUE)
}

if (length(count_files) == 0) {
  stop("No count files found. Please pass a valid count file or run step_2_count.py first.")
}

plot_dir <- file.path("Analysis_Data", "counts", "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

for (file_path in count_files) {
  cat("\n=======================================================\n")
  cat("Processing:", file_path, "\n")

  pool_prefix <- sub("\\.count\\.txt$", "", basename(file_path))
  df <- read.table(file_path, header = TRUE, sep = "\t", check.names = FALSE)

  if (!all(c("sgRNA", "Gene") %in% colnames(df))) {
    warning("File ", file_path, " missing expected sgRNA/Gene columns. Skipping.")
    next
  }

  sample_cols <- setdiff(colnames(df), c("sgRNA", "Gene"))
  # Filter for samples that actually have counts in this library pool
  active_samples <- sample_cols[colSums(df[, sample_cols, drop = FALSE], na.rm = TRUE) > 0]

  if (length(active_samples) == 0) {
    cat("  No active samples with >0 reads in this pool. Skipping.\n")
    next
  }

  cat("  Active samples:", paste(active_samples, collapse = ", "), "\n")

  df_long <- df |>
    pivot_longer(
      cols = all_of(active_samples),
      names_to = "sample",
      values_to = "count"
    )
  
  # Ranked Guide Representation (average count log10 vs ranked gRNA)
  df_ranked <- df |>
    mutate(avg_count = rowMeans(across(all_of(active_samples)), na.rm = TRUE)) |>
    arrange(avg_count) |>
    mutate(rank = row_number())

  p_rep <- ggplot(df_ranked, aes(x = rank, y = log10(avg_count + 1))) +
    geom_area(fill = "purple4", alpha = 0.3) +
    geom_line(color = "purple4", linewidth = 1) +
    labs(
      x = "gRNA Rank (lowest to highest representation)",
      y = "Average Read Count + 1 (log10)",
      title = paste("Guide Representation —", pool_prefix)
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))

  rep_out <- file.path(plot_dir, paste0(pool_prefix, "_guide_representation.png"))
  ggsave(rep_out, p_rep, width = 8, height = 6, dpi = 300)
  cat("  Saved Guide Representation plot:", rep_out, "\n")

  # 1. ECDF of gRNA Representation
  p_ecdf <- ggplot(df_long, aes(x = count + 1, color = sample)) +
    stat_ecdf(geom = "step", linewidth = 1.2) +
    scale_x_log10() +
    labs(
      x = "gRNA Read Count + 1 (log10)",
      y = "Cumulative Fraction of gRNAs",
      title = paste("ECDF of gRNA Representation —", pool_prefix)
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))

  ecdf_out <- file.path(plot_dir, paste0(pool_prefix, "_ecdf.png"))
  ggsave(ecdf_out, p_ecdf, width = 8, height = 6, dpi = 300)
  cat("  Saved ECDF plot:", ecdf_out, "\n")

  # 2. Lorenz Curves + Gini Index
  lorenz_df <- df_long |>
    group_by(sample) |>
    arrange(count) |>
    mutate(
      cum_gRNAs = row_number() / n(),
      cum_reads = cumsum(count) / sum(count)
    ) |>
    ungroup()

  gini_df <- df_long |>
    group_by(sample) |>
    summarize(Gini = round(calc_gini(count), 3), .groups = "drop")

  cat("  Gini Statistics:\n")
  print(gini_df)

  p_lorenz <- ggplot(lorenz_df, aes(x = cum_gRNAs, y = cum_reads, color = sample)) +
    geom_line(linewidth = 1.1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    labs(
      x = "Cumulative Fraction of gRNAs",
      y = "Cumulative Fraction of Total Reads",
      title = paste("Lorenz Curves —", pool_prefix)
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))

  lorenz_out <- file.path(plot_dir, paste0(pool_prefix, "_lorenz.png"))
  ggsave(lorenz_out, p_lorenz, width = 8, height = 6, dpi = 300)
  cat("  Saved Lorenz plot:", lorenz_out, "\n")

  # 3. Dropout Summary Table
  dropout_table <- lorenz_df |>
    group_by(sample) |>
    summarize(
      dropout_1pct = max(count[cum_reads <= 0.01]),
      dropout_5pct = max(count[cum_reads <= 0.05]),
      dropout_10pct = max(count[cum_reads <= 0.10]),
      zero_count_guides = sum(count == 0),
      .groups = "drop"
    )

  cat("  Dropout Summary Table (max read counts at cumulative signal percentiles):\n")
  print(dropout_table)
  dropout_out <- file.path(plot_dir, paste0(pool_prefix, "_dropout_summary.csv"))
  write.csv(dropout_table, dropout_out, row.names = FALSE)

  # 4. Density Overlay
  p_dens <- ggplot(df_long, aes(x = count + 1, fill = sample)) +
    geom_density(alpha = 0.3) +
    scale_x_log10() +
    labs(
      x = "Read count + 1 (log10)",
      y = "Density",
      title = paste("Read Count Density Overlay —", pool_prefix)
    ) +
    theme_bw(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5))

  dens_out <- file.path(plot_dir, paste0(pool_prefix, "_density.png"))
  ggsave(dens_out, p_dens, width = 8, height = 6, dpi = 300)
  cat("  Saved Density plot:", dens_out, "\n")

  # 5. Control Gene gRNA Representation
  ctrl_long <- df_long |>
    filter(grepl("Control|NTC|Non-targeting", Gene, ignore.case = TRUE))

  if (nrow(ctrl_long) > 0 && n_distinct(ctrl_long$Gene) > 0) {
    p_ctrl <- ggplot(ctrl_long, aes(x = count + 1, color = sample)) +
      stat_ecdf(geom = "step", linewidth = 1) +
      scale_x_log10() +
      facet_wrap(~Gene, scales = "free_y") +
      labs(
        x = "Read count + 1 (log10)",
        y = "Cumulative Fraction",
        title = paste("Control gRNA Representation —", pool_prefix)
      ) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5))

    ctrl_out <- file.path(plot_dir, paste0(pool_prefix, "_controls_ecdf.png"))
    ggsave(ctrl_out, p_ctrl, width = 9, height = 7, dpi = 300)
    cat("  Saved Controls plot:", ctrl_out, "\n")
  }
}
cat("\nGuide representation analysis complete.\n")
