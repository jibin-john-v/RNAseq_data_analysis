
# Function to create an MA plot with labeled top genes
create_ma_plot <- function(res, outputfolder,contrast_name) {
    # Create the PDF file
    pdf(glue("{outputfolder}/{contrast_name}_MA_plots.pdf"), height = 12, width = 10)
    # Create the MA plot
    DESeq2::plotMA(res, ylim = c(-5, 5))
    # Add horizontal lines at y = 1 and y = -1 with distinguishable color
    abline(h = 1, col = "blue", lty = 2)  # Horizontal line at y = 1
    abline(h = -1, col = "blue", lty = 2) # Horizontal line at y = -1

    # Sort results by adjusted p-value
    sorted_res <- res[order(res$padj, res$pvalue), ]
    # Identify top 10 upregulated and downregulated genes based on adjusted p-value
    upregulated_genes <- head(rownames(sorted_res[sorted_res$log2FoldChange > 0 & !is.na(sorted_res$padj), ]), 10)
    downregulated_genes <- head(rownames(sorted_res[sorted_res$log2FoldChange < 0 & !is.na(sorted_res$padj), ]), 10)
    # Generate a sequence of offsets for the text labels to avoid overlap
    offset_up <- seq(0.5, 1.5, length.out = length(upregulated_genes))
    offset_down <- seq(-1.5, -0.5, length.out = length(downregulated_genes))

    # Add bold black text labels for top 10 upregulated genes
    for (i in seq_along(upregulated_genes)) {
        gene <- upregulated_genes[i]
        with(sorted_res[gene, ], {
        text(baseMean, log2FoldChange + offset_up[i], gene, col = "black", cex = 0.8, font = 2, adj = c(0.5, 0.5))
        }) }

    # Add bold black text labels for top 10 downregulated genes
    for (i in seq_along(downregulated_genes)) {
        gene <- downregulated_genes[i]
        with(sorted_res[gene, ], {
        text(baseMean, log2FoldChange + offset_down[i], gene, col = "black", cex = 0.8, font = 2, adj = c(0.5, 0.5))
        }) }

    # Close the PDF device
    dev.off()
}



# volcano_plots
deg_enchanced_volcano_plots <- function(res_df2, outputfolder, contrast_name, fc_cutof = 1) {
    # Filter for "Up regulated" and "Down regulated" and sort by FDR
    filtered_genes <- res_df2 %>%
        filter(expression_status_strict %in% c("Up regulated", "Down regulated")) %>%arrange(padj, pvalue)

    # If fewer than 15 genes, supplement with remaining top genes from full res_df2
    if (nrow(filtered_genes) < 15) {
        # Select Up regulated", "Down regulated" based on expression_status
        sorted_genes <- res_df2 %>%
            filter(expression_status %in% c("Up regulated", "Down regulated")) %>%arrange(padj, pvalue)

        # Calculate how many more genes are needed
        remaining_needed <- 15 - nrow(filtered_genes)
        # Take only the remaining genes needed from sorted_genes
        remaining_genes <- head(sorted_genes[!sorted_genes$gene_id %in% filtered_genes$gene_id, ], remaining_needed)
        filtered_genes <- rbind(filtered_genes, remaining_genes)

    }

    # If fewer than 15 genes, supplement with remaining top genes from full res_df2
    if (nrow(filtered_genes) < 15) {
        # Sort the entire dataset by FDR
        sorted_genes <- res_df2 %>% arrange(padj, pvalue)

        # Calculate how many more genes are needed
        remaining_needed <- 15 - nrow(filtered_genes)
        # Take only the remaining genes needed from sorted_genes
        remaining_genes <- head(sorted_genes[!sorted_genes$gene_id %in% filtered_genes$gene_id, ], remaining_needed)
        filtered_genes <- rbind(filtered_genes, remaining_genes)

    }
    # Select top 15 genes for labeling
    selectLab <- head(filtered_genes$gene_id, n = 15)

    # Use tryCatch to handle potential errors without printing warnings or errors
    pvalue_cutoff <- tryCatch(
        {  max_value <- max(res_df2[res_df2$expression_status_strict %in% c("Up regulated", "Down regulated"), ]$pvalue, na.rm = TRUE)
        max_value  # Return the calculated max value
        },error = function(e) {return(NA) },
        warning = function(w) {return(NA)}
    )

    # If pvalue_cutoff is NA, execute the second command
    if (is.na(pvalue_cutoff)) {
        pvalue_cutoff <- res_df2 %>% arrange(padj, pvalue)
        pvalue_cutoff <- min(pvalue_cutoff$pvalue)+0.001
    }

    # Set custom x and y limits
    x_min <- min(res_df2$log2FoldChange) - 1  # Adding buffer
    x_max <- max(res_df2$log2FoldChange) + 1
    y_min <- -log10(max(res_df2$pvalue))  # Use max for PValue to get smallest -log10(PValue)
    y_max <- -log10(min(res_df2$pvalue))  # Use min for PValue to get largest -log10(PValue)

    # Remove existing file if needed
    system(glue('rm {outputfolder}/DGE{contrast_name}_EnhancedVolcano_volcano_plots.pdf') )
    
    # Generate the volcano plot
    pdf(glue("{outputfolder}/{contrast_name}_EnhancedVolcano_volcano_plots.pdf"), height = 8, width = 10)
    print(EnhancedVolcano(res_df2, lab = res_df2$gene_id,
                          x = 'log2FoldChange', y = 'pvalue',
                          pointSize = 2, labSize = 3,
                          selectLab = selectLab,
                          pCutoff = pvalue_cutoff, pCutoffCol = 'pvalue',
                          FCcutoff = fc_cutof, 
                          colAlpha = 0.8,
                          labFace = 'bold',
                          boxedLabels = FALSE,
                          arrowheads=TRUE,
                          drawConnectors = TRUE,
                          widthConnectors = 0.8,
                          lengthConnectors = unit(0.001, "npc"),
                          labCol = 'black',
                          xlim = c(x_min, x_max),  # Adjusted x-axis limits
                          ylim = c(y_min, y_max)   # Custom y-axis limit using -log10(PValue)
    ))
    dev.off()
}


