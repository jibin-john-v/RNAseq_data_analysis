# Get the full script path
path_args <- commandArgs(trailingOnly = FALSE)
#script_dir <- dirname( dirname(normalizePath(sub("--file=", "", path_args[grep("--file=", path_args)]))) )
script_dir <- dirname(normalizePath(sub("--file=", "", path_args[grep("--file=", path_args)])))
source(file.path(script_dir, "visualisation", "deg_visual.R"))
source(file.path(script_dir, "visualisation", "deg_heamap.R"))
source(file.path(script_dir, "visualisation", "exploratory.R"))
source(file.path(script_dir, "visualisation", "heatmap.R"))
source(file.path(script_dir, "visualisation", "pca_mds.R"))
source(file.path(script_dir, "visualisation", "save_figure.R"))



# Function to perform RNA-seq analysis for a specific contrast
dge_separately_contrast_analysis <- function(filtered_countdata_df, selected_metadata_df, filtered_covar_df, filtered_contrast_df, ResultDir, OutputFilename) {
    
    # Ensure the result directory exists
    if (!dir.exists(ResultDir)) {
        dir.create(ResultDir, recursive = TRUE)
    }
    

    if (!is.null(filtered_covar_df)) {
        # Merge metadata and covariates DataFrames
        selected_metadata_df <- merge(
            selected_metadata_df,filtered_covar_df,by = "id")
    
        # Identify all covariate columns (excluding "id" and "condition")
        covariate_columns <- setdiff(colnames(selected_metadata_df), c("id", "condition"))
        # Construct the design formula dynamically
        design_formula <- as.formula(paste("~", paste(c(covariate_columns, "condition"), collapse = " + ")))

    # Create DESeq2 dataset with dynamic covariates
        selected_dds <- DESeqDataSetFromMatrix(
            countData = filtered_countdata_df,
            colData = full_metadata_df[, c("id", "condition", covariate_columns)],
            design = design_formula
        )

    }else {
        selected_metadata_df<-selected_metadata_df
        # Create DESeq2 dataset
        selected_dds <- DESeqDataSetFromMatrix(
        countData = filtered_countdata_df,
        colData = selected_metadata_df[, c("id", "condition")],
        design = ~condition )
    }


    Name <- paste0(filtered_contrast_df$Contrast, "_vs_", filtered_contrast_df$Reference)  # Naming convention
    
    # Estimate size factors and dispersions
    selected_dds <- estimateSizeFactors(selected_dds)
    selected_dds <- estimateDispersions(selected_dds)
    
    # Obtain normalized counts
    selected_normalized_counts <- counts(selected_dds, normalized = TRUE)
    
    # Perform rlog transformation
    selected_rld <- rlog(selected_dds, blind = TRUE)
    
    # Create output directories
    output_directory <- glue("{ResultDir}/DGE/data/{Name}")
    dir.create(output_directory, recursive = TRUE)
    normalized_output_directory <- glue("{ResultDir}/DGE/data/normalized_reads")
    dir.create(normalized_output_directory, recursive = TRUE)
    
    # Save normalized counts
    selected_normalized_countsN <- cbind(Gene_ID = rownames(selected_normalized_counts), selected_normalized_counts)
    write.table(
        selected_normalized_countsN,
        file.path(normalized_output_directory, paste0(OutputFilename, "_ReadCountNormalized.tsv")),
        sep = "\t", row.names = FALSE
    )
    
    # Save rlog-transformed counts
    selected_RlogNCount <- cbind(Gene_ID = rownames(assay(selected_rld)), assay(selected_rld))
    write.table(
        selected_RlogNCount,
        file.path(normalized_output_directory, paste0(OutputFilename, "_ReadCountRlogNorm.tsv")),
        sep = "\t", row.names = FALSE
    )
    
    # Group-wise mean normalized counts
    group_info <- colData(selected_dds)$condition
    selected_mean_groupwise_counts <- data.frame(row.names = rownames(selected_normalized_counts))
    for (group in unique(group_info)) {
        selected_mean_groupwise_counts[[group]] <- rowMeans(selected_normalized_counts[, group_info == group, drop = FALSE])
    }

    # Generate plots and heatmaps
    create_pca_plot(selected_rld, selected_metadata_df, ResultDir, OutputFilename)
    create_mds_plot(selected_rld, selected_metadata_df, ResultDir, OutputFilename)
    
    # Sample-level correlation heatmaps
    create_correlation_heatmap(selected_rld, selected_metadata_df[, c("id", "condition")], ResultDir, glue("sample_level_correlation_", OutputFilename))
    
    # Most variable genes heatmaps
    for (n in c(25, 50, 100, 500, 1000, 2000)) {
        create_most_vargenes_heatmap(n, selected_rld, selected_metadata_df[, c("id", "condition")], ResultDir, OutputFilename)
    }
    
    # Hierarchical clustering dendrogram
    save_rld_hclust_plot_colored(selected_rld, selected_metadata_df, "condition", ResultDir, paste0(OutputFilename, "_Hierarchical_Clustering_Dendrogram"))
    
    # Differential gene expression analysis
    REFERENCE <- as.character(filtered_contrast_df$Reference)
    selected_dds$condition <- relevel(factor(selected_dds$condition), ref = REFERENCE)
    selected_dds <- DESeq(selected_dds)
    
    # Obtain results
    perform_deseq2_analysis(selected_dds, filtered_contrast_df, selected_normalized_counts, selected_metadata_df, selected_mean_groupwise_counts, OutputFilename)
}



# Wrapper function to process all contrasts
separate_deg_all_contrasts <- function(contrast_df, metadata_df, covar_df, countdata_df, ResultDir) {
    for (i in seq_len(nrow(contrast_df))) {
        print(i)  # Print the current contrast index
        CONDITION <- "condition"
        REFERENCE <- toString(contrast_df[i, 1])
        TREATED <- toString(contrast_df[i, 2])  # For naming convention
        OutputFilename <- paste0(TREATED, "_vs_", REFERENCE)
        
        # Filter metadata for the current contrast
        selected_metadata_df <- subset(metadata_df, subset = (condition == TREATED | condition == REFERENCE))
        selected_metadata_df <- selected_metadata_df[, c("id", "condition")]

        # Filter count data for selected samples
        filtered_countdata_df <- countdata_df[, selected_metadata_df$id]

        # Check if column names match IDs (sanity check)
        if (!identical(selected_metadata_df$id, colnames(filtered_countdata_df))) {
        stop("Mismatch between selected_metadata_df IDs and filtered_countdata_df column names!")
        }

        # Filter contrast data for the current pair
        filtered_contrast_df <- as.data.frame(contrast_df[contrast_df$Reference == REFERENCE & contrast_df$Contrast == TREATED, , drop = FALSE])
        
        if (!is.null(covar_df)) {
            filtered_covar_df <- covar_df[covar_df$id %in% selected_metadata_df$id, ]
        }else {
           filtered_covar_df<-NULL
        }

        # Perform RNA-seq analysis for the specific contrast
        dge_separately_contrast_analysis(
        filtered_countdata_df,
        selected_metadata_df,
        filtered_covar_df,
        filtered_contrast_df,
        ResultDir,
        OutputFilename
        )
    }
}


##Deseq analysis function 
perform_deseq2_analysis <- function(dds, contrast_df, normalized_counts, metadata, mean_groupwise_counts,OutputFilename) {
    # Initialize an empty data frame to store merged counts
    all_merged_counts <- data.frame(expression_status = character())
    all_res_df<-data.frame(gene_id = character())

    # Loop over contrasts
    for (i in seq_len(nrow(contrast_df))) {
        
        CONDITION <- "condition"
        REFERENCE <- toString(contrast_df[i, 1])
        TREATED <- toString(contrast_df[i, 2])
        Name <- paste0(TREATED,"_vs_",REFERENCE)
        
        # DESeq2 results
        res <- results(dds, contrast = c(CONDITION, TREATED, REFERENCE))
        
        selected_metadata <- subset(metadata, subset = (condition == TREATED | condition == REFERENCE))
        selected_samples <- selected_metadata$id
        selected_normalized_counts <- subset(normalized_counts, select = selected_samples)
        
        res_df <- as.data.frame(res)
        res_df$expression_status_strict <- "No change"
        res_df$expression_status_strict[res$log2FoldChange >= 1 & res$padj < 0.05] <- "Up regulated"
        res_df$expression_status_strict[res$log2FoldChange <= -1 & res$padj < 0.05] <- "Down regulated"
        res_df$expression_status_strict[res$log2FoldChange >= 1 & res$padj >= 0.05] <- "Not significant"
        res_df$expression_status_strict[res$log2FoldChange <= -1 & res$padj >= 0.05] <- "Not significant"
        
        res_df$expression_status <- "No change"
        res_df$expression_status[res$log2FoldChange >= 0 & res$padj < 0.05] <- "Up regulated"
        res_df$expression_status[res$log2FoldChange <= 0 & res$padj < 0.05] <- "Down regulated"
        res_df$expression_status[res$log2FoldChange >= 0 & res$padj >= 0.05] <- "Not significant"
        res_df$expression_status[res$log2FoldChange <= 0 & res$padj >= 0.05] <- "Not significant"
        res_df <- res_df[order(res_df$padj), ]

        res_df_data_NormReadCount <- merge(res_df, mean_groupwise_counts, by = 0, all = TRUE)
        colnames(res_df)[colnames(res_df) == 'Row.names'] <- 'gene_id'
        res_df$gene_id <- row.names(res_df)

        colnames(res_df_data_NormReadCount)[colnames(res_df_data_NormReadCount) == 'Row.names'] <- 'gene_id'
        
        # Create directory for results
        #dir.create(paste0(glue("DES/DGE/data/{Name}"), Name), recursive = TRUE, showWarnings = FALSE)
        write.csv(res_df_data_NormReadCount, file.path(paste0(ResultDir,"/DGE/data/", Name), paste0(OutputFilename, "_DESeq2_with_NormReadCount.csv")), row.names = FALSE) 
        write.csv(res_df, file.path(paste0(ResultDir,"/DGE/data/", Name), paste0(OutputFilename, "_DESeq2.csv")), row.names = FALSE) 
        
        # Count unique values for expression_status_strict
        merged_counts <- full_join(
            res_df %>%
                dplyr::count(expression_status_strict) %>%
                rename(expression_status = expression_status_strict, !!glue("{Name}_Strict") := n),
            res_df %>%dplyr::count(expression_status) %>%
                rename(!!glue("{Name}_Lenient") := n), 
            by = "expression_status" )
        
        # Store merged counts
        all_merged_counts <- full_join(all_merged_counts, merged_counts, by = "expression_status")
        # Create plots
        dir.create(file.path(glue('{ResultDir}/DGE/plots/{Name}') ), recursive = TRUE, showWarnings = FALSE)
        deg_enchanced_volcano_plots(res_df, glue('{ResultDir}/DGE/plots/{Name}'), OutputFilename)
        create_ma_plot(res, glue('{ResultDir}/DGE/plots/{Name}'), OutputFilename)
        # Heat maps
        deg_heatmap(res_df, 50, normalized_counts, selected_metadata[, c('id', 'condition')], glue('{ResultDir}/DGE/plots/{Name}'), OutputFilename)
        deg_heatmap(res_df, 100, normalized_counts, selected_metadata[, c('id', 'condition')], glue('{ResultDir}/DGE/plots/{Name}'), OutputFilename)
        deg_heatmap(res_df, 100, normalized_counts, selected_metadata[, c('id', 'condition')], glue('{ResultDir}/DGE/plots/{Name}'), OutputFilename)
        deg_heatmap(res_df, 100, normalized_counts, selected_metadata[, c('id', 'condition')], glue('{ResultDir}/DGE/plots/{Name}'), OutputFilename)
        deg_heatmap(res_df, 100, normalized_counts, selected_metadata[, c('id', 'condition')], glue('{ResultDir}/DGE/plots/{Name}'), OutputFilename)

        # Add prefix to all columns except gene_id
        res_df <- res_df %>%rename_with(~ ifelse(. == "gene_id", ., paste(OutputFilename, ., sep = "_")), -gene_id)
        all_res_df <- full_join(all_res_df, res_df, by = "gene_id")
    }

    # Save all merged counts to a single CSV file
    dir.create(glue("{ResultDir}/DGE/data/all_conditions_merged_results"))
    write.csv(all_merged_counts, file.path(glue("{ResultDir}/DGE/data/all_conditions_merged_results"), glue("{OutputFilename}_DGE_Genes_counts.csv")  ), row.names = FALSE)
    write.csv(all_res_df, file.path(glue("{ResultDir}/DGE/data/all_conditions_merged_results"), glue("{OutputFilename}_DGE_Genes_merged.csv") ), row.names = FALSE)

    strict_dge_df <- all_res_df %>%select(gene_id, matches("expression_status_strict"))
    lenient_dge_df <- all_res_df %>%select(gene_id, matches("expression_status"), -matches("_status_strict$"))
    write.csv(strict_dge_df, file.path(glue("{ResultDir}/DGE/data/all_conditions_merged_results"), glue("{OutputFilename}_DGE_Strict_Genes_merged_SelectedColumns.csv") ), row.names = FALSE)
    write.csv(lenient_dge_df, file.path(glue("{ResultDir}/DGE/data/all_conditions_merged_results"), glue("{OutputFilename}_DGE_Lenient_Genes_merged_SelectedColumns.csv") ), row.names = FALSE)

    return(all_res_df)
}


##Interaction Deseq to analysis 
perform_deseq2_interaction_analysis <- function(int_dds,countdata_df, metadata_df) {
    ###Preparing separate normalised read for savin to fle and 
    int_normalized_counts <- counts(int_dds, normalized=TRUE)
    # calculate groupwise (condition) mean
    group_info <- colData(int_dds)$condition
    mean_groupwise_counts <- data.frame(row.names = rownames(int_normalized_counts))
    # Calculate mean normalized counts for each group
    for (group in unique(group_info)) {
        mean_groupwise_counts[[group]] <- rowMeans(int_normalized_counts[, group_info == group, drop = FALSE])
        }# Create a data frame to hold the results

    ### Transform counts for data visualization: https:vscode-file://vscode-app/private/var/folders/rl/w43l6j1j5_s67l6k8gn5hlqm0000gp/T/AppTranslocation/54C764F7-C9F7-40BA-816B-6E688787A1DD/d/VisualStudioCode.app/Contents/Resources/app/out/vs/code/electron-sandbox/workbench/workbench.html//hbctraining.github.io/DGE_workshop_salmon/lessons/03_DGE_QC_analysis.html; 
            #Transform normalized counts using the rlog transformation
    #vsd <- vst(dds, blind=FALSE)
    int_rld <- rlog(int_dds, blind=TRUE)  ## assay() is function from the "SummarizedExperiment" package that was loaded when you loaded DESeq2
    interaction_term<-resultsNames(int_dds)[length(resultsNames(int_dds))]

    # Initialize an empty data frame to store merged counts and results
    all_merged_counts <- data.frame(expression_status = character())
    all_res_df <- data.frame(gene_id = character())
    Name <- paste0("Interaction_", interaction_term)

    # DESeq2 results
    res <- results(int_dds, name = interaction_term)
    res_df <- as.data.frame(res)

    # Add expression status columns
    res_df$expression_status_strict <- "No change"
    res_df$expression_status_strict[res$log2FoldChange >= 1 & res$padj < 0.05] <- "Up regulated"
    res_df$expression_status_strict[res$log2FoldChange <= -1 & res$padj < 0.05] <- "Down regulated"
    res_df$expression_status_strict[res$log2FoldChange >= 1 & res$padj >= 0.05] <- "Not significant"
    res_df$expression_status_strict[res$log2FoldChange <= -1 & res$padj >= 0.05] <- "Not significant"
    
    res_df$expression_status <- "No change"
    res_df$expression_status[res$log2FoldChange >= 0 & res$padj < 0.05] <- "Up regulated"
    res_df$expression_status[res$log2FoldChange <= 0 & res$padj < 0.05] <- "Down regulated"
    res_df$expression_status[res$log2FoldChange >= 0 & res$padj >= 0.05] <- "Not significant"
    res_df$expression_status[res$log2FoldChange <= 0 & res$padj >= 0.05] <- "Not significant"

    # Merge with normalized group-wise counts
    res_df_data_NormReadCount <- merge(res_df, mean_groupwise_counts, by = "row.names", all = TRUE)
    colnames(res_df_data_NormReadCount)[colnames(res_df_data_NormReadCount) == 'Row.names'] <- 'gene_id'

    # Ensure gene_id is the first column in res_df
    res_df$gene_id <- rownames(res_df)
    res_df <- res_df[, c("gene_id", setdiff(colnames(res_df), "gene_id"))]

    # Create directory for results
    dir.create(paste0(ResultDir,"/DGE/data/", Name), recursive = TRUE, showWarnings = FALSE)

    # Save results to CSV
    write.csv(res_df_data_NormReadCount, file.path(paste0(ResultDir,"/DGE/data/", Name), paste0(Name, "_DESeq2_with_NormReadCount.csv")), row.names = FALSE) 
    write.csv(res_df, file.path(paste0(ResultDir,"/DGE/data/", Name), paste0(Name, "_DESeq2.csv")), row.names = FALSE) 
    
    # Count unique values for expression_status_strict and expression_status
    merged_counts <- full_join(
        res_df %>%count(expression_status_strict) %>%
            rename(expression_status = expression_status_strict, !!glue("{Name}_Strict") := n),
        res_df %>%
            count(expression_status) %>%rename(!!glue("{Name}_Lenient") := n), by = "expression_status" )

    # Create plots directory
    dir.create(file.path(glue(ResultDir,'/DGE/plots/{Name}') ), recursive = TRUE, showWarnings = FALSE)
    
    # Generate volcano and MA plots
    deg_enchanced_volcano_plots(res_df_data_NormReadCount, glue(ResultDir,'/DGE/plots/{Name}'), Name)
    create_ma_plot(res, glue(ResultDir,'/DGE/plots/{Name}'), Name)

    # Generate heatmaps with different gene counts
    for (n_genes in c(50, 100,500)) {
        deg_heatmap(res_df_data_NormReadCount, n_genes, int_normalized_counts, metadata_df[, c('id', 'condition')], glue(ResultDir,'/DGE/plots/{Name}'), Name)
    }

    # Add prefix to all columns except gene_id
    res_df <- res_df %>% rename_with(~ ifelse(. == "gene_id", ., paste(Name, ., sep = "_")), -gene_id)
    all_res_df <- full_join(all_res_df, res_df, by = "gene_id")

    # Save all merged counts and results
    write.csv(merged_counts, file.path(glue('{ResultDir}/DGE/data/all_conditions_merged_results/'), glue("{interaction_term}_Interaction_DGE_Genes_counts.csv") ), row.names = FALSE)
    write.csv(all_res_df, file.path(glue('{ResultDir}/DGE/data/all_conditions_merged_results/'), glue("{interaction_term}_Interaction_DGE_Genes_merged.csv") ), row.names = FALSE)

    # Save strict and lenient DGE summaries
    strict_dge_df <- all_res_df %>% select(gene_id, matches("expression_status_strict"))
    lenient_dge_df <- all_res_df %>% select(gene_id, matches("expression_status"), -matches("_status_strict$"))
    
    write.csv(strict_dge_df, file.path(glue('{ResultDir}/DGE/data/all_conditions_merged_results/'), glue("{interaction_term}_Interaction_DGE_Strict_Genes_merged_SelectedColumns.csv"  ) ), row.names = FALSE)
    write.csv(lenient_dge_df, file.path(glue('{ResultDir}/DGE/data/all_conditions_merged_results/'), glue("{interaction_term}_Interaction_DGE_Lenient_Genes_merged_SelectedColumns.csv" ) ), row.names = FALSE)

    return(all_res_df)
}
