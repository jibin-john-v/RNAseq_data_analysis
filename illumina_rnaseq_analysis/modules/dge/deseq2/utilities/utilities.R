
## Process count data 
prepare_count_data <- function(countdata) {
  
  # Convert input to a data frame if it's not already
  countdata_df <- as.data.frame(countdata)
  
  # Get the gene ID column name
  geneid_column <- colnames(countdata_df)[1]
  
  # Set row names using the gene ID column
  row.names(countdata_df) <- countdata_df[[geneid_column]]
  
  # Remove the gene ID column from the data frame
  countdata_df <- countdata_df[, -1]
  
  minimum_count<-nrow(countdata_df)/2
  # Filter out rows (genes) with total read counts below the minimum count
  countdata_df <- countdata_df[rowSums(countdata_df) >= minimum_count, ]
  
  return(countdata_df)
}



## This function will check both covar_df and metadata_df contain id column
check_id_column <- function(covar_df, metadata_df) {
    # Initialize an empty vector to collect errors
    error_messages <- c()
    
    # Check if 'id' column exists in covar_df
    if (!is.null(covar_df)) {
            if (!"id" %in% colnames(covar_df)) {
        error_messages <- c(error_messages, "Error: The 'id' column is missing in the covariate data frame (covar_df).")
        }
    } 
    
    # Check if 'id' column exists in metadata_df
    if (!"id" %in% colnames(metadata_df)) {
        error_messages <- c(error_messages, "Error: The 'id' column is missing in the metadata data frame (metadata_df).")
    }

    # If any errors were found, terminate the program and display all error messages
    if (length(error_messages) > 0) {
        stop(paste(error_messages, collapse = "\n"))
    }

}


# Function to create the 'condition' column
prepare_meta_data <- function(metadata_df, id_col='id') {
    # Identify columns to concatenate (all columns after the `id_col`)
    id_index <- which(colnames(metadata_df) == id_col)
    if (length(id_index) == 0) {
      stop("The specified ID column does not exist in the data frame.")
    }
    # Concatenate all columns after the ID column
    metadata_df <- metadata_df %>%
      mutate(condition = apply(df[, (id_index + 1):ncol(df)], 1, function(row) paste(row, collapse = "_")))
        
  return(metadata_df)
}


##Check input count data sample name and metadata sample name are in same order or not
check_inputdata <- function(countdata_df, metadata_df,covar_df) {
    metadata_df <- metadata_df %>%mutate_all(~ gsub("-", "_", .))
    countdata_df <- countdata_df %>%mutate_all(~ gsub("-", "_", .))
    countdata_df <- `rownames<-`(data.frame(lapply(countdata_df, as.numeric)), rownames(countdata_df))

    # Check if the column names in countdata match the metadata IDs
    if (!is.null(covar_df)) { 
        if (!all(sort(colnames(countdata_df)) == sort(metadata_df$id)) || 
            !all(sort(colnames(countdata_df)) == sort(covar_df$id)) || 
            !all(sort(metadata_df$id) == sort(covar_df$id))) {
            stop("Error: Column names do not match between countdata, metadata, and covar data.")
        }
    } else {
        if (!all(sort(colnames(countdata_df)) == sort(metadata_df$id))) {
            stop("Error: Column names do not match between countdata and metadata.")
        }
    }
    
    # Sort metadata
    metadata_df <- metadata_df[order(metadata_df$condition, metadata_df$id), ]
    # Reorder countdata columns
    countdata <- countdata[, metadata$id]

    if (!is.null(covar_df)) {
        covar_df <- covar_df[match(metadata_df$id, covar_df$id), ]
        }

    if (!is.null(covar_df)) {
        # Final check to ensure alignment
        if (!all(colnames(countdata_df) == metadata_df$id & metadata_df$id == covar_df$id)) {
            stop("\n\nSample names in the ReadCount file, SampleInfo file and covar file are not the same or not in the same order; please correct..\n\n")
        }
    }else{
        if (!all(colnames(countdata_df) == metadata_df$id )) {
            stop("\n\nSample names in the ReadCount file and SampleInfo files are not the same or not in the same order; please correct..\n\n")
        }
    }
    return(list(countdata_df = countdata_df, metadata_df = metadata_df, covar_df=covar_df))
}


contrast_check <- function(metadata_df,contrast_df) {
    # Check if all values in contrast_df$Reference are present in metadata_df$condition
    if (!all(contrast_df$Reference %in% metadata_df$condition)) {
        stop("Error: Some values in contrast_df$Reference are not present in metadata_df$condition.")
    } 
    if (!all(contrast_df$Contrast %in% metadata_df$condition)) {
        stop("Error: Some values in contrast_df$Contrast are not present in metadata_df$condition.")
    } 
return(contrast_df)
}