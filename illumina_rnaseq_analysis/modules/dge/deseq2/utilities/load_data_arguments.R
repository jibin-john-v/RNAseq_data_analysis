# Define a function to parse arguments and load files
load_data_and_arguments <- function() {
    # Argument parsing
    parser <- ArgumentParser()
    parser$add_argument("--read_count",
        required = TRUE,
        help = "Provide read count file, it should only contain feature name and counts")
        
    parser$add_argument("--sample_info", 
        required = TRUE,
        help = "Provide file with sample info; headers [id]")

    parser$add_argument("--covar_file", 
        required = FALSE,  # Make this argument optional
        default = NULL,    # Default to NULL if not provided
        help = "Provide file with sample covariates file; headers [id,covariates] (optional)")

    parser$add_argument("--contrast_file", 
        required = TRUE,
        help = "Provide file with contrast info; headers [Reference,Contrast]")

    parser$add_argument("--output", 
        required = TRUE,
        help = "Provide Output folder name")
    
    # Parse command-line arguments
    args <- parser$parse_args()
    
    # Load files using provided arguments
    countdata_df <- fread(args$read_count)
    metadata_df <- fread(args$sample_info)
    contrast_df <- fread(args$contrast_file)
    
    # Load covariate file only if provided
    if (!is.null(args$covar_file)) {
        covar_df <- fread(args$covar_file)
    } else {
        covar_df <- NULL  # Set to NULL if not provided
    }
    
    OutputFilename <- args$output
    
    # Return a list of loaded data
    return(list(
        countdata_df = countdata_df,
        metadata_df = metadata_df,
        contrast_df = contrast_df,
        covar_df = covar_df,  # Can be NULL if not provided
        OutputFilename = OutputFilename
    ))
}
