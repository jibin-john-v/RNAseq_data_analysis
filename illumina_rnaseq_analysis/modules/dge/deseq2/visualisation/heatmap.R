
# Define a function to create heatmaps and save in all formats
create_correlation_heatmap <- function(rld, metadata, ResultDir, OutputFilename,out_dir, 
                                       cluster_rows = TRUE, cluster_columns = TRUE, 
                                       row_split = NULL, column_split = NULL, row_gap = 0.2, column_gap = 0.2) {
  
    dir.create(file.path(ResultDir, "summary","sample_correlation",out_dir), recursive = TRUE, showWarnings = FALSE)
    out_path <- glue("{ResultDir}/summary/sample_correlation/",out_dir) 

    # Extract the rlog matrix from the rld object
    rld_mat <- assay(rld)
    # Compute pairwise correlation values
    rld_cor <- cor(rld_mat)
    
    # Prepare the annotation colors
        annotation_col <- select(metadata, condition)
        unique_conditions <- unique(annotation_col$condition)
        condition_colors <- viridis::viridis(length(unique_conditions))
        annotation_colors <- setNames(condition_colors, unique_conditions)

    # Create heatmap annotations for rows and columns
    ha_col <- HeatmapAnnotation(condition = annotation_col$condition,
                                col = list(condition = annotation_colors), show_legend = FALSE)
    ha_row <- rowAnnotation(condition = annotation_col$condition,
                            col = list(condition = annotation_colors), show_legend = TRUE)
    
    # Define a list of file types and their corresponding graphic device functions
    file_types <- c("tiff", "png", "pdf", "svg")
    device_functions <- list(
        tiff = function() tiff(file = file.path(out_path, paste0(OutputFilename, ".tiff")), units = "in", width = 8, height = 8, res = 800),
        png = function() png(file = file.path(out_path, paste0(OutputFilename, ".png")), units = "in", width = 8, height = 8, res = 800),
        pdf = function() pdf(file = file.path(out_path, paste0(OutputFilename, ".pdf")), width = 8, height = 8),
        svg = function() svg(file = file.path(out_path, paste0(OutputFilename, ".svg")), width = 8, height = 8)
    )
    
    # Loop over each file type and save the heatmap in all formats
    for (file_type in file_types) {
        # Open the corresponding graphic device
        device_functions[[file_type]]()
        # Create the heatmap object
        ht <- Heatmap(rld_cor,
                    top_annotation = ha_col,left_annotation = ha_row,
                    cluster_rows = cluster_rows, cluster_columns = cluster_columns,   # Cluster columns based on user input
                    row_split = row_split, column_split = annotation_col$condition,   # Optionally split columns
                    row_gap = unit(row_gap, "mm"),   # Row gap
                    column_gap = unit(column_gap, "mm"),   # Column gap
                    rect_gp = gpar(col = "white", lwd = 0.1),  # Set border color and line width
                    show_heatmap_legend = TRUE,   # Show heatmap legend
                    column_dend_height = unit(2, "cm"), 
                    row_dend_width = unit(2, "cm"),
                    heatmap_legend_param = list(direction = "vertical", 
                                                title = "correlation", title_position = "topcenter", 
                                                title_gp = gpar(fontsize = 10, fontface = "bold")))
        # Draw the heatmap
        draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
        
        # Close the graphic device
        dev.off()
    }
}



# Function to save hierarchical clustering plot for DESeq2 rld object with color-coded dendrogram
save_rld_hclust_plot_colored <- function(rld, metadata, condition_column = "condition", output_folder, output_filename, file_types = c("pdf", "png", "tiff", "svg"), width = 8, height = 8, res = 300) {
    # Create the output directory for saving results
    dir.create(file.path(output_folder, "summary",'sample_hclust',output_filename), recursive = TRUE, showWarnings = FALSE)
    out_path <- file.path(output_folder, "summary",'sample_hclust',output_filename)
    
    # Extract assay data from the rld object and perform hierarchical clustering
    assay_data <- assay(rld)

    # Create distance matrix and hierarchical clustering object
    hclust_obj <- hclust(dist(t(assay_data), method = "euclidean"))

    # Create a dendrogram object
    dend <- as.dendrogram(hclust_obj)

    # Check that metadata contains the condition column and matches colnames of assay data
    if (!condition_column %in% colnames(metadata)) {
        stop(paste("Condition column", condition_column, "not found in metadata."))
    }
    if (!all(colnames(assay_data) %in% metadata$id)) {
        stop("Sample names in assay data and metadata do not match.")
    }

    # Order metadata to match the columns of assay_data
    metadata <- metadata[match(colnames(assay_data), metadata$id), ]

    # Assign colors to samples based on the condition column
    conditions <- as.factor(metadata[[condition_column]])
    n_conditions <- length(levels(conditions))

    # Generate a color palette that avoids fluorescent colors
    if (n_conditions <= 8) {
        # Use the "Dark2" palette for up to 8 conditions (avoids fluorescent colors)
        condition_colors <- brewer.pal(n_conditions, "Dark2")
    } else {
        # For more conditions, use a custom set of colors
        condition_colors <- c(
            "#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", 
            "#A6761D", "#666666", "#8B0000", "#4682B4", "#6B8E23", "#DA70D6",
            "#B22222", "#00008B", "#9932CC", "#556B2F", "#8FBC8F"
        )
        # Recycle colors if conditions exceed the length of the palette
        condition_colors <- rep(condition_colors, length.out = n_conditions)
    }
    condition_colors <- setNames(condition_colors, levels(conditions))

    # Label the branches of the dendrogram with colors based on condition
    labels_colors(dend) <- condition_colors[conditions][order.dendrogram(dend)]

    # Define device functions for each file type
    device_functions <- list(
        pdf = function(file) pdf(file = file, width = width, height = height),
        png = function(file) png(file = file, width = width, height = height, units = "in", res = res),
        tiff = function(file) tiff(file = file, width = width, height = height, units = "in", res = res),
        svg = function(file) svg(file = file, width = width, height = height)
    )

    # Loop through each specified file type and create the corresponding plot
    for (file_type in file_types) {
        if (file_type %in% names(device_functions)) {
            # Construct the output file name
            output_file <- file.path(out_path, paste0(output_filename, ".", file_type))

            # Open the appropriate graphics device
            device_functions[[file_type]](output_file)

            # Create the hierarchical clustering plot with colored branches
            plot(dend, main = "Hierarchical Clustering Dendrogram", ylab = "Height")

            # Add legend for conditions
            legend("topright", legend = levels(conditions), col = condition_colors, pch = 19, title = condition_column)

            # Close the graphics device
            dev.off()
        } else {
            warning(paste("File type", file_type, "is not supported. Skipping..."))
        }
    }
}

