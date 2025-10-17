analyze_tf_enrichment <- function(ipa_gr, polya_gr, tf_dir, output_dir, 
                                  exp_name, bkd_name, threshold_for_TFs_Enrichment,
                                  min_overlap = 1) {
  cat("=== Starting TF Enrichment Analysis ===\n\n")
  tf_files <- list.files(tf_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(tf_files) == 0) {
    stop("No .rds files found in the specified directory")
  }
  cat(sprintf("Found %d TF files to process\n", length(tf_files)))
  cat(sprintf("%s ranges: %d\n", exp_name, length(ipa_gr)))
  cat(sprintf("%s ranges: %d\n", bkd_name, length(polya_gr)))
  cat(sprintf("Distance threshold: %s\n", threshold_for_TFs_Enrichment))
  cat("\n")
  n_ipa <- length(ipa_gr)
  n_polya <- length(polya_gr)
  results_list <- list()
  for (i in seq_along(tf_files)) {
    tf_file <- tf_files[i]
    tf_name <- gsub("\\.rds$", "", basename(tf_file))
    cat(sprintf("[%d/%d] Processing TF: %s\n", i, length(tf_files), tf_name))
    tryCatch({
      tf_gr <- readRDS(tf_file)
      if (!is(tf_gr, "GRanges")) {
        warning(sprintf("  Skipping %s: not a GRanges object", tf_name))
        next
      }
      cat(sprintf("  TF binding sites: %d\n", length(tf_gr)))
      overlaps_ipa <- findOverlaps(ipa_gr, tf_gr, minoverlap = min_overlap)
      observed_ipa <- length(unique(queryHits(overlaps_ipa)))
      overlaps_polya <- findOverlaps(polya_gr, tf_gr, minoverlap = min_overlap)
      observed_polya <- length(unique(queryHits(overlaps_polya)))
      cat(sprintf("  %s with TF binding: %d\n", exp_name, observed_ipa))
      cat(sprintf("  %s with TF binding: %d\n", bkd_name, observed_polya))
      p_background <- observed_polya / n_polya
      expected_ipa <- n_ipa * p_background
      variance_ipa <- n_ipa * p_background * (1 - p_background)
      if (variance_ipa == 0) {
        z_score <- NA
        p_value <- NA
        cat("  Warning: variance is zero, skipping z-score calculation\n")
      } else {
        sd_ipa <- sqrt(variance_ipa)
        z_score <- (observed_ipa - expected_ipa) / sd_ipa
        p_value <- pnorm(z_score, lower.tail = FALSE)
        cat(sprintf("  Expected in %s: %.2f\n", exp_name, expected_ipa))
        cat(sprintf("  Z-score: %.4f\n", z_score))
        cat(sprintf("  P-value: %.4e\n", p_value))
      }
      results_list[[tf_name]] <- data.frame(
        TF_name = tf_name,
        Observed_in_Exp = observed_ipa,
        Observed_in_Bkd = observed_polya,
        Expected_in_Exp = expected_ipa,
        Background_probability = p_background,
        Z_score = z_score,
        P_value = p_value,
        stringsAsFactors = FALSE
      )
      cat("\n")
    }, error = function(e) {
      warning(sprintf("  Error processing %s: %s", tf_name, e$message))
      cat("\n")
    })
  }
  if (length(results_list) == 0) {
    stop("No TFs were successfully processed")
  }
  results_df <- bind_rows(results_list)
  results_df <- results_df %>%
    arrange(P_value)
  results_df$P_value_adjusted <- p.adjust(results_df$P_value, method = "BH")
  results_df$Significant <- ifelse(results_df$P_value_adjusted < 0.05, "Yes", "No")
  cat("=== Analysis Complete ===\n")
  cat(sprintf("Total TFs analyzed: %d\n", nrow(results_df)))
  cat(sprintf("Significant TFs (FDR < 0.05): %d\n", 
              sum(results_df$Significant == "Yes", na.rm = TRUE)))
  output_file <- file.path(output_dir, sprintf("TF_enrichment_%s_%s_%s.csv", 
                                               exp_name, bkd_name, threshold_for_TFs_Enrichment))
  write.csv(results_df, output_file, row.names = FALSE)
  cat(sprintf("\nResults saved to: %s\n", output_file))
  output_rds <- file.path(output_dir, sprintf("TF_enrichment_%s_%s_%s.rds", 
                                              exp_name, bkd_name, threshold_for_TFs_Enrichment))
  saveRDS(results_df, output_rds)
  cat(sprintf("Results also saved as RDS: %s\n", output_rds))
  return(results_df)
}

run_tf_enrichment_analysis <- function(ipa_gr, polya_gr, tf_dir, output_dir, 
                                       exp_name = "IPA", bkd_name = "PolyA",
                                       threshold_for_TFs_Enrichment = "500bp",
                                       significance_threshold = 0.05) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  results <- analyze_tf_enrichment(ipa_gr, polya_gr, tf_dir, output_dir, 
                                   exp_name, bkd_name, threshold_for_TFs_Enrichment)
  results$Significant <- ifelse(results$P_value_adjusted < significance_threshold, "Yes", "No")
  cat(sprintf("\nUsing significance threshold: %.4f\n", significance_threshold))
  cat(sprintf("Significant TFs at this threshold: %d\n", 
              sum(results$Significant == "Yes", na.rm = TRUE)))
  results_plot <- results %>%
    filter(!is.na(Z_score) & !is.na(P_value))
  if (nrow(results_plot) > 0) {
    p_volcano <- ggplot(results_plot, aes(x = Z_score, y = -log10(P_value))) +
      geom_point(aes(color = Significant), alpha = 0.6, size = 2) +
      geom_hline(yintercept = -log10(significance_threshold), 
                 linetype = "dashed", color = "red", size = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", size = 0.8) +
      scale_color_manual(values = c("Yes" = "#F8766D", "No" = "gray60")) +
      labs(
        title = sprintf("TF Enrichment in %s vs %s with %s distance around", 
                        exp_name, bkd_name, threshold_for_TFs_Enrichment),
        x = "Z-score",
        y = "-log10(P-value)",
        color = paste0("Significant\n(FDR < ", significance_threshold, ")")
      ) +
      theme_classic() +
      theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 11, color = "black"),
        axis.title = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.background = element_blank(),
        legend.title = element_text(face = "bold"),
        axis.ticks = element_line(color = "black")
      )
    volcano_file <- file.path(output_dir, sprintf("TF_enrichment_%s_%s_%s_volcano_plot.pdf", 
                                                  exp_name, bkd_name, threshold_for_TFs_Enrichment))
    pdf(volcano_file, width = 8, height = 8)
    print(p_volcano)
    dev.off()
    cat(sprintf("\nVolcano plot saved to: TF_enrichment_%s_%s_%s_volcano_plot.pdf\n", 
                exp_name, bkd_name, threshold_for_TFs_Enrichment))
  }
  top_n <- min(20, nrow(results_plot))
  top_tfs <- results_plot %>%
    arrange(P_value) %>%
    head(top_n)
  if (nrow(top_tfs) > 0) {
    p_bar <- ggplot(top_tfs, aes(x = reorder(TF_name, Z_score), y = Z_score, fill = Significant)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_manual(values = c("Yes" = "#F8766D", "No" = "gray60")) +
      labs(
        title = sprintf("Top %d TFs by Enrichment (%s vs %s with %s distance around)", 
                        top_n, exp_name, bkd_name, threshold_for_TFs_Enrichment),
        x = "Transcription Factor",
        y = "Z-score",
        fill = paste0("Significant\n(FDR < ", significance_threshold, ")")
      ) +
      theme_classic() +
      theme(
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.background = element_blank(),
        axis.line = element_line(color = "black", size = 0.5),
        axis.text = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 9),
        axis.title = element_text(size = 12, color = "black", face = "bold"),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "right",
        legend.background = element_blank(),
        legend.title = element_text(face = "bold"),
        axis.ticks = element_line(color = "black")
      )
    bar_file <- file.path(output_dir, sprintf("TF_enrichment_%s_%s_%s_top_barplot.pdf", 
                                              exp_name, bkd_name, threshold_for_TFs_Enrichment))
    pdf(bar_file, width = 8, height = 8)
    print(p_bar)
    dev.off()
    cat(sprintf("Bar plot saved to: TF_enrichment_%s_%s_%s_top_barplot.pdf\n", 
                exp_name, bkd_name, threshold_for_TFs_Enrichment))
  }
  output_file <- file.path(output_dir, sprintf("TF_enrichment_%s_%s_%s.csv", 
                                               exp_name, bkd_name, threshold_for_TFs_Enrichment))
  write.csv(results, output_file, row.names = FALSE)
  output_rds <- file.path(output_dir, sprintf("TF_enrichment_%s_%s_%s.rds", 
                                              exp_name, bkd_name, threshold_for_TFs_Enrichment))
  saveRDS(results, output_rds)
  return(results)
}