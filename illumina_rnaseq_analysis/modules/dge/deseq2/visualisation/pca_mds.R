

# PCA Plotting Function
create_pca_plot <- function(deseq_object, selected_metadata_df,ResultDir,OutputFilename) {
    # Generate ellipse data
    draw_ellipse <- function(df) {
        group_means <- aggregate(cbind(PC1, PC2) ~ condition, data = df, mean)
        ellipses <- list()
        for (i in 1:nrow(group_means)) {
            group_data <- df[df$condition == group_means$condition[i], ]
            cov_matrix <- cov(group_data[, c("PC1", "PC2")])
            ellipses[[i]] <- as.data.frame(ellipse(cov_matrix, centre = as.numeric(group_means[i, 2:3]), level = 0.95))
            ellipses[[i]]$condition <- group_means$condition[i]
        }
        do.call(rbind, ellipses)
    }
    dir.create(file.path(ResultDir, "summary",'pca',{OutputFilename}), recursive = TRUE, showWarnings = FALSE)
    out_path <- glue("{ResultDir}/summary/pca/{OutputFilename}/")
    
    # Create DESeq2 dataset and rlog transform
    pcaData <- plotPCA(deseq_object, intgroup = c("condition","id"), returnData = TRUE, ntop = 1000)
    percentVar <- round(100 * attr(pcaData, "percentVar"))
    covariates <- setdiff(colnames(selected_metadata_df), c("id", "condition"))

    if (length(covariates) > 0) {
        for (covariate in covariates) {
            toplot <- data.frame(PC1 = pcaData$PC1, PC2 = pcaData$PC2, condition = selected_metadata_df$condition, covariate = selected_metadata_df[[covariate]])
            ellipse_data <- draw_ellipse(toplot)
            p <- ggplot(toplot, aes(PC1, PC2, colour = condition, shape = covariate)) +
                geom_point(size = 1.5) +
                geom_polygon(data = ellipse_data, aes(x = PC1, y = PC2, fill = condition), alpha = 0.25) +
                labs(x = glue("PC1 ({percentVar[1]}%)"),y = glue("PC2 ({percentVar[2]}%)"),
                     color = "Condition",shape = covariate) +theme_minimal() +
                theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
                      axis.text = element_text(size = 12),
                      axis.title = element_text(size = 14, face = "bold"))
            # Save the current plot
            ggsave(glue("{out_path}{OutputFilename}_PCA_{covariate}_PC1_PC2.pdf"), plot = p, device = "pdf", width = 12, height = 8)
            print(glue("PCA Plot with covariate '{covariate}' saved."))
        }
    } else {
        toplot <- data.frame(PC1 = pcaData$PC1, PC2 = pcaData$PC2, condition = selected_metadata_df$condition, sample_name = selected_metadata_df$id)
        ellipse_data <- draw_ellipse(toplot)
        p <- ggplot(toplot, aes(PC1, PC2, colour = condition, shape = condition)) +geom_point(size = 1) +
            geom_polygon(data = ellipse_data, aes(x = PC1, y = PC2, fill = condition), alpha = 0.25) +
            labs(x = glue("PC1 ({percentVar[1]}%)"),y = glue("PC2 ({percentVar[2]}%)")) +theme_minimal() +
            theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.text = element_text(size = 12),axis.title = element_text(size = 14, face = "bold"))
        save_plot(p, out_path, glue('{OutputFilename}_PCA1_2_with_ellipse'), width = 12, height = 8, dpi = 300)

        p2 <- ggplot(toplot, aes(PC1, PC2, colour = condition, shape = condition)) +
            geom_point(size = 2) +
            labs(x = glue("PC1 ({percentVar[1]}%)"),y = glue("PC2 ({percentVar[2]}%)")) +theme_minimal() +
            theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.text = element_text(size = 12),axis.title = element_text(size = 14, face = "bold"))
        save_plot(p2, out_path, glue('{OutputFilename}_PCA1_2_without_ellipse_without_label'), width = 12, height = 8, dpi = 300)


        p3 <- ggplot(toplot, aes(PC1, PC2, colour = condition, shape = condition)) +
            geom_point(size = 2) +
            geom_text_repel(aes(label = sample_name), 
                            size = 3, max.overlaps = Inf, # Allow all labels (or set a specific number)
                            box.padding = 0.5,       # Space around the text box
                            point.padding = 0.5,     # Space around the points
                            nudge_x = 0.1,           # Slight nudge on x-axis
                            nudge_y = 0.1,           # Slight nudge on y-axis
                            segment.color = 'grey50' # Color for the line connecting the point and label
                            ) + 
            labs(x = glue("PC1 ({percentVar[1]}%)"), y = glue("PC2 ({percentVar[2]}%)")) +
            theme_minimal() +theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
            axis.text = element_text(size = 12),axis.title = element_text(size = 14, face = "bold"))
        save_plot(p3, out_path, glue('{OutputFilename}_PCA1_2_without_ellipse_with_label'), width = 12, height = 8, dpi = 300)
    }
}


# MDS Plotting Function
create_mds_plot <- function(selected_rld, selected_metadata_df, outputfolder,OutputFilename) {
    out_path <- glue("{ResultDir}/summary/mds/{OutputFilename}/") 
    dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
    expression_data <- assay(selected_rld)

    # Generate MDS plot data
    mds <- limma::plotMDS(expression_data)
    covariates <- setdiff(colnames(selected_metadata_df), c("id", "condition"))

    if (length(covariates) > 0) {
        for (covariate in covariates) {
            toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, condition = selected_metadata_df$condition, covariate = selected_metadata_df[[covariate]])

            p <- ggplot(toplot, aes(Dim1, Dim2, colour = condition, shape = covariate)) +
                geom_point(size = 3) +
                labs(x = glue("MDS1 ({round(mds$var.explained[1], 2) * 100}%)"),
                     y = glue("MDS2 ({round(mds$var.explained[2], 2) * 100}%)"),
                     color = "Condition",shape = covariate) + theme_minimal()

            # Save the current plot
            ggsave(glue("{out_path}{OutputFilename}_MDS_{covariate}.pdf"), plot = p, device = "pdf", width = 8, height = 6)
            print(glue("MDS Plot with covariate '{covariate}' saved."))
        }
    } else {
        toplot <- data.frame(Dim1 = mds$x, Dim2 = mds$y, condition = selected_metadata_df$condition, sample_name = selected_metadata_df$id)
        p1 <- ggplot(toplot, aes(Dim1, Dim2, colour = condition,shape =condition )) +
            geom_point(size = 3) +
            labs(x = glue("MDS1 ({round(mds$var.explained[1], 2) * 100}%)"),
                 y = glue("MDS2 ({round(mds$var.explained[2], 2) * 100}%)"),color = "Condition") +
            theme_minimal()
        save_plot(p1, out_path, glue('{OutputFilename}_MDS_Plot_no_label'), width = 12, height = 8, dpi = 300)

        # Create the plot with sample names added
        p2 <- ggplot(toplot, aes(Dim1, Dim2, colour = condition, shape = condition)) +
            geom_point(size = 3) +
            geom_text_repel(aes(label = sample_name), size = 3, max.overlaps = Inf, box.padding = 0.5, point.padding = 0.5) +  # Add sample names
            labs(x = glue("MDS1 ({round(mds$var.explained[1], 2) * 100}%)"),
                y = glue("MDS2 ({round(mds$var.explained[2], 2) * 100}%)"), 
                color = "Condition") + theme_minimal()
            save_plot(p2, out_path, glue('{OutputFilename}_MDS_Plot_with_label'), width = 12, height = 8, dpi = 300)

    }
}
