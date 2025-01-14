# Rscript src/tapsformer/differential_methylation/dss_prepare_rastair_data.r --suffix raw
# Rscript src/tapsformer/differential_methylation/dss_prepare_rastair_data.r --suffix raw_with_liver

start_time <- Sys.time()

install_and_load <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    install.packages(new.pkg, dependencies = TRUE)
  }
  sapply(pkg, require, character.only = TRUE)
}

# Function for Bioconductor packages
bioc_install_and_load <- function(pkg) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) {
    BiocManager::install(new.pkg, update = FALSE)
  }
  sapply(pkg, require, character.only = TRUE)
}

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Define packages
bioc_packages <- c(
  "DSS", "GenomicRanges", "bsseq", "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "AnnotationHub"
)
cran_packages <- c("data.table", "futile.logger", "parallel", "dplyr", "optparse")

# Install and load packages
bioc_install_and_load(bioc_packages)
install_and_load(cran_packages)

# Print loaded package versions
sessionInfo()

# Parse command line arguments
option_list <- list(
  make_option(c("-s", "--suffix"),
    type = "character", default = "raw",
    help = "Suffix for base directory [default= %default]", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Construct base_dir with the provided suffix
base_dir <- file.path("/users/zetzioni/sharedscratch/tapsformer/data/methylation/by_cpg", opt$suffix)

# Set up logging
flog.appender(appender.file("/users/zetzioni/sharedscratch/logs/dss_prepare.log"))

flog.info("Finished loading libraries for DSS analysis script")
flog.info("Starting DSS analysis script")
flog.info(paste("Using base directory:", base_dir))
# Set up logging
flog.appender(appender.file("/users/zetzioni/sharedscratch/logs/dss_prepare.log"))

flog.info("Finished loading libraries for DSS analysis script")
flog.info("Starting DSS analysis script")

num_cores <- 12
# Create a cluster
cl <- makeCluster(num_cores)

get_sample_name <- function(file_path, prefix) {
  # Extract the base file name without the path
  sample_name <- basename(file_path)

  # Remove the '_ScrBsl_rastair.bed' or '_Ctrl_rastair.bed' suffix to extract the <name>
  sample_name <- gsub("_ScrBsl_rastair\\.bed$|_Ctrl_rastair\\.bed$", "", sample_name)

  # Create the final RDS file name with the prefix
  final_name <- paste0(prefix, "_", sample_name, ".rds")

  return(final_name)
}

# Function to read and preprocess bed files - we are aggregating the cpgs from the output of rastair where they are split
# across the C>T and G>A
read_preprocess_bed <- function(file_path) {
  flog.info(paste("Reading and preprocessing file:", file_path))
  if (!file.exists(file_path)) {
    flog.error(paste("File does not exist:", file_path))
    return(NULL)
  }
  dt <- fread(file_path, header = TRUE, skip = 0)
  if (nrow(dt) == 0) {
    flog.warn(paste("File is empty:", file_path))
    return(NULL)
  }
  flog.info(paste("Original columns:", paste(names(dt), collapse = ", ")))
  flog.info(paste("Original number of rows:", nrow(dt)))
  if ("#chr" %in% names(dt)) {
    setnames(dt, "#chr", "chr")
  }
  # Filter for chromosomes 1-22
  dt <- dt[chr %in% paste0("chr", 1:22)]
  flog.info(paste("Rows after filtering for chr1-22:", nrow(dt)))
  # Create a unique identifier for each CpG site
  dt[, cpg_id := paste(chr, ifelse(strand == "+", end, start), sep = "_")]
  # Count unique CpG sites before merging
  unique_cpgs <- uniqueN(dt$cpg_id)
  flog.info(paste("Number of unique CpG sites:", unique_cpgs))
  # Ensure position represents the start of the CpG
  dt[, pos := ifelse(strand == "+", start, end)] # Start for + strand, end for - strand
  # Combine data for each CpG site
  result <- dt[, .(
    chr = chr[1],
    pos = min(pos), # Use the minimum position to represent the CpG start
    N = sum(unmod + mod),
    X = sum(mod) # X represents methylated cytosines (mod in TAPS)
  ), by = .(cpg_id)]

  # Calculate beta
  result[, `:=`(
    beta = X / N, # beta represents proportion of methylated cytosines
    cpg_id = NULL
  )]
  if (nrow(result) == 0) {
    flog.warn(paste("No data left after filtering for file:", file_path))
    return(NULL)
  }
  flog.info(paste("Finished preprocessing", file_path))
  flog.info(paste("Number of unique CpG sites after preprocessing:", nrow(result)))
  flog.info(paste("Columns in preprocessed data:", paste(names(result), collapse = ", ")))

  return(result)
}

# Read tumour and control files
flog.info("Reading tumour and control files")
flog.info(paste("Using base directory:", base_dir))
flog.info("tumour file pattern: ScrBsl.*\\.bed$")
flog.info("Control file pattern: Ctrl.*\\.bed$")

tumour_files <- list.files(path = base_dir, pattern = "ScrBsl.*\\.bed$", full.names = TRUE) # nolint
control_files <- list.files(path = base_dir, pattern = "Ctrl.*\\.bed$", full.names = TRUE) # nolint

flog.info(paste("Found", length(tumour_files), "tumour files and", length(control_files), "control files")) # nolint

# Log the file names for verification
flog.info("tumour files:")
sapply(tumour_files, function(f) flog.info(paste("  ", f)))
flog.info("Control files:")
sapply(control_files, function(f) flog.info(paste("  ", f)))

# Define a wrapper function that includes logging
read_preprocess_bed_wrapper <- function(file_path) {
  library(data.table)
  library(futile.logger)

  # Set up logging for each worker
  flog.appender(appender.file("/users/zetzioni/sharedscratch/logs/dss_analysis_worker.log"), name = "worker")

  result <- tryCatch(
    {
      read_preprocess_bed(file_path)
    },
    error = function(e) {
      flog.error(paste("Error processing file:", file_path, "-", conditionMessage(e)), name = "worker")
      return(NULL)
    }
  )

  return(result)
}

# Export the necessary functions and objects to the cluster
clusterExport(cl, c("read_preprocess_bed_wrapper", "read_preprocess_bed"))

# Process tumour files in parallel
save_sample_data <- function(data_list, base_dir, prefix, file_paths) {
  for (i in seq_along(data_list)) {
    sample_data <- data_list[[i]]
    if (!is.null(sample_data)) {
      # Get the sample name and format it for saving
      sample_rds_name <- get_sample_name(file_paths[i], prefix)
      sample_rds_path <- file.path(base_dir, sample_rds_name)

      # Save the data
      saveRDS(sample_data, sample_rds_path)
      flog.info(paste("Sample", sample_rds_name, "data saved to:", sample_rds_path))
    }
  }
}

# Process tumour files in parallel
flog.info("Processing tumour files in parallel")
tumour_data_list <- parLapply(cl, tumour_files, read_preprocess_bed_wrapper)

# Process control files in parallel
flog.info("Processing control files in parallel")
control_data_list <- parLapply(cl, control_files, read_preprocess_bed_wrapper)

# Stop the cluster
stopCluster(cl)

# Remove NULL elements (if any) from the results
tumour_data_list <- tumour_data_list[!sapply(tumour_data_list, is.null)]
control_data_list <- control_data_list[!sapply(control_data_list, is.null)]

# Save the data as separate RDS files for each sample
flog.info("Saving individual tumour samples to RDS")
save_sample_data(tumour_data_list, base_dir, "tumour", tumour_files)

flog.info("Saving individual control samples to RDS")
save_sample_data(control_data_list, base_dir, "control", control_files)

flog.info("Preparation script completed successfully")

end_time <- Sys.time()
flog.info(paste("Total runtime:", difftime(end_time, start_time, units = "mins"), "minutes"))
