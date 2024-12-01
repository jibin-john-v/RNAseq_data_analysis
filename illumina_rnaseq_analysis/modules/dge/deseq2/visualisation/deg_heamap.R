

##Create DEG Heat map 
deg_heatmap <- function(res_df, n_genes, normalized_counts, metadata, outputfolder,contrast_name) {

        # Order metadata by condition 
        metadata<-metadata[order(metadata$condition), ]
        
        # Filter for "Up regulated" and "Down regulated" and sort by FDR and p-value
        filtered_df <- res_df %>%
            filter(expression_status_strict %in% c("Up regulated", "Down regulated")) %>%arrange(padj, pvalue)

        # If fewer than n_genes are found, supplement with remaining top genes
        if (nrow(filtered_df) < n_genes) {
            # Select Up and Down regulated genes based on expression_status
            sorted_genes <- res_df %>%filter(expression_status %in% c("Up regulated", "Down regulated")) %>% arrange(padj, pvalue)
            # Calculate how many more genes are needed
            remaining_needed <- n_genes - nrow(filtered_df)
            # Select remaining genes from sorted_genes that are not already in filtered_df
            remaining_genes <- head(sorted_genes[!sorted_genes$gene_id %in% filtered_df$gene_id, ], remaining_needed)
            # Combine filtered_df with remaining_genes, ensuring there are no duplicates
            filtered_df <- rbind(filtered_df, remaining_genes)
        }

        # If filtered_df is still empty, or fewer than n_genes, take top genes by FDR from full res_df
        if (nrow(filtered_df) == 0) {
            message("No 'Up regulated' or 'Down regulated' genes found. Selecting top genes by FDR.")
            filtered_df <- res_df %>% arrange(padj, pvalue) %>% head(n = n_genes)
        } else if (nrow(filtered_df) < n_genes) {
            # Supplement with top genes from res_df if fewer than n_genes are found
            remaining_genes <- res_df[!res_df$gene_id %in% filtered_df$gene_id, ] %>%
                arrange(padj, pvalue) %>%head(n = n_genes - nrow(filtered_df))
            filtered_df <- rbind(filtered_df, remaining_genes)
        }

        # Finally, take only the top n_genes if there are more than needed
        filtered_df <- head(filtered_df, n = n_genes)
        # Check for duplicates
        filtered_df <- filtered_df[!duplicated(filtered_df$gene_id), ]

    # Define device functions for each file type
    device_functions <- list(
        tiff = function(file) tiff(file = file, units = "in", width = 8, height = 8, res = 800),
        png = function(file) png(file = file, units = "in", width = 8, height = 8, res = 800),
        pdf = function(file) pdf(file = file, width = 8, height = 8),
        svg = function(file) svg(file = file, width = 8, height = 8) )
    
    file_types = c("tiff", "png", "pdf", "svg")
    # Select rows from log2_cpm_counts based on filtered_df
    selected_rows <- normalized_counts[rownames(normalized_counts) %in% filtered_df$gene_id, ]
    selected_rows <- selected_rows[, metadata$id]

    metadata2 <- metadata[order(metadata$condition), ]
    selected_rows <- selected_rows[, metadata2$id]
    
    # Generate color palettes for annotations
    conditions <- unique(metadata2$condition)
    n_conditions <- length(conditions)
    
    # Use ggsci or RColorBrewer for distinct colors
    condition_colors <- pal_npg("nrc")(n_conditions)  # Use a palette from ggsci
    names(condition_colors) <- conditions

    # Prepare annotation colors for other metadata columns if needed
    ann_colors <- list(condition = condition_colors)

    # Create the annotation object with differentiable colors
    #ann <- data.frame(subset(metadata2, select = -id))
    #colnames(ann) <- colnames(metadata2)[-1]
    #colAnn <- HeatmapAnnotation(df = ann, col = ann_colors)

        # Prepare annotation colors
        annotation_col <- metadata[, "condition", drop = FALSE]
        unique_conditions <- unique(annotation_col$condition)
        condition_colors <- viridis::viridis(length(unique_conditions))
        annotation_colors <- setNames(condition_colors, unique_conditions)

        # Create column and row annotations
        ha_col <- HeatmapAnnotation(
            condition = annotation_col$condition,
            col = list(condition = annotation_colors), show_legend = TRUE
        )

        ha_row <- rowAnnotation(
            condition = annotation_col$condition,
            col = list(condition = annotation_colors), show_legend = FALSE
        )

        ##Added 1e-6 insted of 0, to avoid error
        selected_rows[selected_rows == 0] <- 1e-6
    
        # Create the heatmap object
        ht <- Heatmap(
            t(scale(t(selected_rows))),
            name = "z-score",
            top_annotation = ha_col,
            column_split = metadata2$condition,
            show_row_names = n_genes <= 50,
            show_column_names = nrow(metadata) <= 50,
            show_column_dend=TRUE,
            column_dend_reorder = TRUE,
            rect_gp = gpar(col = "white", lwd = 0), #if not want gap then coment it
            border = FALSE,
            row_names_gp = gpar(fontsize = 8),
            show_heatmap_legend = TRUE,  # Enable single legend for heatmap
            heatmap_legend_param = list(
                direction = "vertical",
                title = "z-score",
                title_position = "topcenter",
                title_gp = gpar(fontsize = 10, fontface = "bold")
            ),
            show_row_dend = n_genes <= 50,
            row_gap = unit(10, "mm"),
        )

        # Loop through each file type and create the corresponding output
        for (file_type in file_types) {
            if (file_type %in% names(device_functions)) {
                # Construct the file name
                file_name <- file.path(outputfolder, paste0(contrast_name,"_Top_", n_genes, "_Genes_heatmap.", file_type))
                # Open the appropriate graphics device
                device_functions[[file_type]](file_name)

                # Draw the heatmap with legends merged into a single column
                draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right", merge_legend = TRUE)

                # Close the graphics device
                dev.off()
            } else {
                warning(paste("File type", file_type, "is not supported. Skipping..."))
            }
        }
    }