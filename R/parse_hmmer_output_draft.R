library(dplyr)
library(readr)
library(stringr)
library(fs)

# Set paths
species <- "Efa"
input_dir <- "hmmer_results"
mapping_file <- file.path("hmmer_inputs", paste0(species, "_feature_sequence_mapping.tsv"))
output_file <- file.path("hmmer_results", paste0(species, "_annotated_features.tsv"))

# Load mapping table
mapping <- read_tsv(mapping_file, show_col_types = FALSE)

# Function to parse HMMer tabular output
parse_hmmer_tbl <- function(tbl_file, source) {
  lines <- read_lines(tbl_file)
  lines <- lines[!str_starts(lines, "#")]
  if (length(lines) == 0) {
    return(tibble())
  }

  df <- read_table2(paste(lines, collapse = "\n"), col_names = FALSE, comment = "#")
  colnames(df)[1:4] <- c("target_name", "target_accession", "query_name", "query_accession")
  df <- df %>%
    group_by(query_name) %>%
    slice_min(order_by = X5, n = 1, with_ties = FALSE) %>% # X5 is E-value
    ungroup() %>%
    mutate(source = source)
  return(df)
}

# Parse all .tbl files
tbl_files <- dir_ls(input_dir, regexp = "\\.tbl$")
parsed_hits <- bind_rows(lapply(tbl_files, function(f) {
  source <- str_remove(path_file(f), "\\.tbl$")
  parse_hmmer_tbl(f, source)
}))

# Join with mapping table
annotated <- mapping %>%
  rename(query_name = Variable) %>%
  left_join(parsed_hits, by = "query_name")

# Save final annotated table
write_tsv(annotated, output_file)
message("Annotated feature table saved as ", output_file)
