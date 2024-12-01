
# Function to generate enhanced boxplot, density plots, and histograms
generate_sample_summary_plots <- function(normalized_counts, metadata, ResultDir, OutputFilename) {

    dir.create(file.path(ResultDir, "summary","sample_levelQC",OutputFilename), recursive = TRUE, showWarnings = FALSE)
    out_path <- glue("{ResultDir}/summary/sample_levelQC/",OutputFilename) 

    OutputFilename=glue("{OutputFilename}_sample_level_QC")

    # List of file formats to save plots
    formats <- c("tiff", "pdf", "png", "svg")
        
    # Function to save a plot in all specified formats
    save_plot_in_formats <- function(plot_function, file_prefix) {
        for (file_format in formats) {
        file_path <- file.path(out_path, paste0(file_prefix, ".", file_format))
        # Open the graphics device for the specified format
        switch(file_format,
                tiff = tiff(file_path, units = "in", width = 12, height = 12, res = 800),
                pdf = pdf(file_path, width = 12, height = 12),
                png = png(file_path, units = "in", width = 12, height = 12, res = 800),
                svg = svg(file_path, width = 12, height = 12))
        # Call the provided plotting function
        plot_function()
        dev.off()
        }
    }

    # Enhanced Boxplot with sorting by 75th percentile
    boxplot_function <- function() {
        # Calculate the 75th percentile (3rd quartile) for each sample (column) of log2(normalized counts)
        percentile_25 <- apply(log2(normalized_counts), 2, quantile, probs = 0.25)
        sorted_indices <- order(percentile_25)
        sorted_counts <- log2(normalized_counts)[, sorted_indices]
        # Define color palette for the sorted samples
        box_colors <- colorRampPalette(c("lightblue", "blue", "darkblue"))(ncol(sorted_counts))
        # Create the boxplot
        boxplot(sorted_counts, 
                col = box_colors,          # Add color to the boxes
                border = "darkblue",       # Darker border color for contrast
                outline = FALSE,           # Hide outliers to reduce clutter
                notch = TRUE,              # Add notches for visual comparison
                pch = 19,                  # Shape of points for outliers
                horizontal = TRUE,         # Horizontal layout
                cex.axis = 0.8,            # Increase axis font size for better readability
                las = 1,                   # Rotate labels for better alignment
                boxwex = 0.7,              # Box width
                ylab = "Samples",          # Label for the Y-axis
                xlab = "log2(normalized Counts)",  # Label for the X-axis
                main = "Distribution of Gene Counts Across Samples (Sorted by 75th Percentile)",  # Add a main title
                cex.main = 1.2,            # Increase the title size
                cex.lab = 1.5,             # Increase the size of axis labels
                font.lab = 2)              # Increase the title size
        # Add grid lines for better clarity
        grid(nx = NULL, ny = NULL, col = "gray90", lty = "dotted")
    }

    # Density Plots
    density_plot_function <- function() {
        plotDensity(log2(normalized_counts), 
                    xlab = "log2(normalized counts)", 
                    cex.lab = 0.7, 
                    panel.first = grid()) 
        legend("topright", legend = c(metadata$id), lwd = 2)
    }
    
    # Histograms of counts per gene
    histogram_function <- function() {
        hist(as.matrix(log2(normalized_counts)), 
            breaks = 100, 
            col = "blue", 
            border = "white",
            main = "Log2-transformed normalized counts per gene", 
            xlab = "log2(counts+1)", 
            ylab = "Number of genes", 
            las = 1, 
            cex.axis = 0.7)
    }

    # Save each type of plot in all formats
    save_plot_in_formats(boxplot_function, paste0(OutputFilename, "_Gene_Count_DistributionsBoxplots"))
    save_plot_in_formats(density_plot_function, paste0(OutputFilename, "_Density"))
    save_plot_in_formats(histogram_function, paste0(OutputFilename, "_Count_PerGeneHist"))
}





# Define a function to create heatmaps and save in all formats
# Function to create correlation heatmaps and save in multiple formats
create_correlation_heatmap <- function(rld, metadata, result_dir, output_filename, out_dir, 
                                       cluster_rows = TRUE, cluster_columns = TRUE, 
                                       row_split = NULL, column_split = NULL, 
                                       row_gap = 0.2, column_gap = 0.2) {
  
    # Prepare output directory
    heatmap_dir <- file.path(result_dir, "summary", "sample_correlation", out_dir)
    dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Path for saving heatmaps
    out_path <- file.path(heatmap_dir, output_filename)

    # Extract the rlog matrix and compute pairwise correlation
    rld_mat <- assay(rld)
    rld_cor <- cor(rld_mat)


    # Prepare annotations
    annotation_col <- select(metadata, condition)
    annotation_colors <- generate_annotation_colors(annotation_col$condition)
    

    ha_col <- HeatmapAnnotation(
        condition = annotation_col$condition,
        col = list(condition = annotation_colors),
        show_legend = FALSE
    )
    
    ha_row <- rowAnnotation(
        condition = annotation_col$condition,
        col = list(condition = annotation_colors),
        show_legend = TRUE
    )
    
    # Define graphic devices for file formats
    graphic_devices <- list(
        tiff = function(file) tiff(file, units = "in", width = 8, height = 8, res = 800),
        png = function(file) png(file, units = "in", width = 8, height = 8, res = 800),
        pdf = function(file) pdf(file, width = 8, height = 8),
        svg = function(file) svg(file, width = 8, height = 8)
    )
    
    # Generate and save heatmaps
    save_heatmap <- function(suffix, row_split, column_split) {
        for (file_type in names(graphic_devices)) {
        file_name <- paste0(out_path, suffix, ".", file_type)
        graphic_devices[[file_type]](file_name)
        
        # Create and draw heatmap
        ht <- Heatmap(
            rld_cor,
            top_annotation = ha_col,
            left_annotation = ha_row,
            cluster_rows = cluster_rows,
            cluster_columns = cluster_columns,
            row_split = row_split,
            column_split = column_split,
            row_gap = unit(row_gap, "mm"),
            column_gap = unit(column_gap, "mm"),
            rect_gp = gpar(col = "white", lwd = 0.1),
            show_heatmap_legend = TRUE,
            column_dend_height = unit(2, "cm"),
            row_dend_width = unit(2, "cm"),
            heatmap_legend_param = list(
            direction = "vertical",
            title = "correlation",
            title_position = "topcenter",
            title_gp = gpar(fontsize = 10, fontface = "bold")
            )
        )
        draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
        
        # Close the graphic device
        dev.off()
        }
    }
    
    # Save heatmap with splits
    save_heatmap("_split", row_split, annotation_col$condition)
    
    # Save heatmap without splits
    save_heatmap("_nosplit", NULL, NULL)
    }

    # Helper function to generate annotation colors
    generate_annotation_colors <- function(conditions) {
    unique_conditions <- unique(conditions)
    condition_colors <- viridis::viridis(length(unique_conditions))
    setNames(condition_colors, unique_conditions)
}




# Define the updated heatmap function
create_most_vargenes_heatmap <- function(n_genes, rld, metadata, outputfolder = outputName,output_filename, file_types = c("tiff", "png", "pdf", "svg")) {

    metadata<-metadata[order(metadata$condition), ]
    metadata2 <- metadata[order(metadata$condition), ]

    # Create the output directory for saving results
    dir.create(file.path(outputfolder, "summary",'var_genes_heat_map',output_filename), recursive = TRUE, showWarnings = FALSE)
    out_path <- file.path(outputfolder, "summary",'var_genes_heat_map',output_filename)
    
    # Define device functions for each file type
    device_functions <- list(
        tiff = function(file) tiff(file = file, units = "in", width = 8, height = 8, res = 800),
        png = function(file) png(file = file, units = "in", width = 8, height = 8, res = 800),
        pdf = function(file) pdf(file = file, width = 8, height = 8),
        svg = function(file) svg(file = file, width = 8, height = 8)
    )

    # Extract assay data from the rld object
    input_matrix <- assay(rld)
    
    # Calculate variance for each gene and select the top N variable genes
    var_genes <- apply(input_matrix, 1, var)
    select_var <- names(sort(var_genes, decreasing = TRUE))[1:n_genes]
    highly_variable_lcpm <- input_matrix[select_var, ]

    # Order metadata and subset based on condition
    col_split <- metadata$condition
    highly_variable_lcpm <- highly_variable_lcpm[, metadata$id]

    # Check that sample names match between assay data and metadata
    if (any(colnames(highly_variable_lcpm) != metadata$id)) {
        stop("\n\nSample names in the ReadCount file and SampleInfo file are not the same or not in the same order; please correct..\n\n")
    } else {
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


    if ( n_genes <500 ){
        # Create the heatmap object
        ht <- Heatmap(
            t(scale(t(highly_variable_lcpm))),
            name = "z-score",
            top_annotation = ha_col,
            column_split = col_split,  # Split columns based on condition
            cluster_rows = n_genes <= 500,
            cluster_columns = TRUE,
            show_column_dend=FALSE,
            show_row_names = n_genes <= 50,
            show_column_names = nrow(metadata) <= 50,
            rect_gp = gpar(col = "white", lwd = 0.001),
            border = FALSE,
            show_heatmap_legend = TRUE,  # Enable single legend for heatmap
            heatmap_legend_param = list(
                direction = "vertical",title = "z-score",
                title_position = "topcenter",
                title_gp = gpar(fontsize = 10, fontface = "bold")),show_row_dend = n_genes <= 50 )

    }else {
            ht <- Heatmap(
            t(scale(t(highly_variable_lcpm))),
            name = "z-score",
            top_annotation = ha_col,
            column_split = col_split,  # Split columns based on condition
            cluster_rows =TRUE,cluster_columns = TRUE,
            show_column_dend=TRUE,
            show_row_names = n_genes <= 50,
            show_column_names = nrow(metadata) <= 50,
            border = FALSE,
            show_heatmap_legend = TRUE,  # Enable single legend for heatmap
            heatmap_legend_param = list(
                direction = "vertical",title = "z-score",
                title_position = "topcenter",
                title_gp = gpar(fontsize = 10, fontface = "bold")),show_row_dend = n_genes <= 50)
    }

        # Loop through each file type and create the corresponding output
        for (file_type in file_types) {
            if (file_type %in% names(device_functions)) {
                # Construct the file name
                file_name <- file.path(out_path, paste0(output_filename,"_Top_", n_genes, "_Variable_Genes_heatmap.", file_type))

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
