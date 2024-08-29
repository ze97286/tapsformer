start_time <- Sys.time()

source("dss_common.r")
plot_dmls_per_dmr <- function(filtered_dmls_in_dmrs, combined_bsseq, output_dir, top_n_dmrs = 20) {
  # Get unique DMRs and sort by areaStat to get top N
  top_dmrs <- unique(filtered_dmls_in_dmrs[order(-abs(as.numeric(dmr_areaStat)))][1:top_n_dmrs, .(dmr_id, dmr_areaStat)])
  
  # Get the raw methylation data
  methylation_data <- bsseq::getMeth(combined_bsseq, type = "raw")
  sample_names <- colnames(methylation_data)
  
  for (i in 1:nrow(top_dmrs)) {
    current_dmr <- top_dmrs[i]
    dmr_dmls <- filtered_dmls_in_dmrs[dmr_id == current_dmr$dmr_id]
    
    plot_data <- data.table(chr = dmr_dmls$chr, pos = dmr_dmls$start)
    
    for (j in 1:nrow(dmr_dmls)) {
      dml <- dmr_dmls[j]
      matching_indices <- which(seqnames(combined_bsseq) == dml$chr & start(combined_bsseq) == dml$start)
      meth_levels <- methylation_data[matching_indices, , drop = FALSE]
      
      if (!is.null(meth_levels) && nrow(meth_levels) > 0 && ncol(meth_levels) > 0) {
        plot_data[j, (sample_names) := as.list(meth_levels)]
      }
    }
    
    # Melt the data for plotting
    plot_data <- melt(plot_data,
      id.vars = c("chr", "pos"),
      variable.name = "Sample", value.name = "MethylationLevel"
    )
    
    plot_data[, SampleType := ifelse(grepl("^tumour", Sample), "Tumour", "Control")]
    
    if (nrow(plot_data) == 0) {
      print(sprintf("No plot data for DMR %s", current_dmr$dmr_id))
      next
    }
    
    # Generate and save the plot
    output_filename <- file.path(output_dir, sprintf("dml_methylation_plot_dmr_%s.svg", gsub(":", "_", current_dmr$dmr_id)))
    
    tryCatch({
      dir.create(dirname(output_filename), showWarnings = FALSE, recursive = TRUE)
      
      num_loci <- nrow(dmr_dmls)
      ncol <- min(5, num_loci)
      nrow <- ceiling(num_loci / ncol)
      
      plot_width <- max(18, ncol * 2)
      plot_height <- max(10, nrow * 1.5)
      
      svglite::svglite(output_filename, width = plot_width, height = plot_height)
      
      plot <- ggplot(plot_data, aes(x = Sample, y = MethylationLevel, shape = SampleType, color = SampleType)) +
        geom_point(size = 2) +
        facet_wrap(~ chr + pos, scales = "free_y", ncol = ncol) +
        theme_bw() +
        labs(
          title = sprintf("Methylation Levels for DMLs in DMR %s (areaStat: %s)", current_dmr$dmr_id, current_dmr$dmr_areaStat),
          x = "Sample", y = "Methylation Level"
        ) +
        scale_shape_manual(values = c(Tumour = 16, Control = 1)) +
        scale_color_manual(values = c(Tumour = "blue", Control = "red")) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
          strip.text = element_text(size = 12),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          panel.spacing = unit(1.5, "lines")
        )
      
      print(plot)
      dev.off()
      
      print(sprintf("Plot saved for DMR %s", current_dmr$dmr_id))
    },
    error = function(e) {
      print(sprintf("Error generating plot for DMR %s: %s", current_dmr$dmr_id, conditionMessage(e)))
      if (!is.null(dev.list())) dev.off()
    })
  }
}

main <- function() {
    # Set paths
    base_dir <- "/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg/subset3/dmr_delta_0.50_p_0.0100_fdr_0.01_minCpG_4_minLen_50_disMerge_50_smooth_cv_window"
    output_dir <- file.path(base_dir, "dml_plots")
    combined_bsseq_file <- file.path(base_dir, "combined_bsseq.rds")
    print(combined_bsseq_file)
    dml_file <- file.path(base_dir, "dmls_in_high_confidence_hypomethylated_dmrs.bed")
    print(dml_file)

    # Load combined BSseq object
    combined_bsseq <- readRDS(combined_bsseq_file)

    # Load DML BED file
    dml_data <- fread(dml_file,
        header = TRUE,
        col.names = c(
            "chr", "start", "end", "name", "score", "strand",
            "dml_pval", "dml_stat", "dmr_areaStat", "dmr_id", "strength"
        )
    )
    setDT(dml_data)

    # Call the plotting function with dml_data instead of filtered_dmls_in_dmrs
    plot_dmls_per_dmr(dml_data, combined_bsseq, output_dir, top_n_dmrs = 20)
}
# Run the main function
main()
