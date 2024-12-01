set.seed("1234")
suppressMessages(library(argparse))

# Get the full script path
path_args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(normalizePath(sub("--file=", "", path_args[grep("--file=", path_args)])))
# Construct the path to utilities.R and source it


source(file.path(script_dir, "utilities", "utilities.R"))
source(file.path(script_dir, "utilities", "load_data_arguments.R"))
source(file.path(script_dir, "deseq_dge", "dge.R"))

config_params<-yaml::read_yaml(file.path(script_dir, "deseq2_config.yaml") )



##Load all the 
data <- load_data_and_arguments()

# Access data from the returned list
OutputFolder <- data$OutputFilename



##check and Prepare the input data
countdata_df<-prepare_count_data(data$countdata_df) 

check_id_column(data$covar_df, data$metadata_df) 

metadata_df<-prepare_meta_data(metadata_df, id_col='id')

inputdata_dfs=check_inputdata(countdata_df, metadata_df,data$covar_df)
countdata_df<-inputdata_dfs$countdata_df
metadata_df<-inputdata_dfs$metadata_df
covar_df<-inputdata_dfs$covar_df

contrast_df<-contrast_check(metadata_df,data$contrast_df)


separate_deg_all_contrasts(contrast_df, metadata_df, countdata_df, ResultDir)