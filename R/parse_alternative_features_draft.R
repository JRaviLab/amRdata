### Testing export to Parquet
library(tidyverse)
library(arrow)

### --- Set me! ---

# Which bug are we doin'?
bug <- "Aba"

# Which feature are we doin'?
feature <- "ARG" # COG or ARG mapping file

# Running this from the repo or from Alpine?
alpine <- TRUE

### --- This stuff should run unsupervised ---

# Determine the input and output path stuff based on the variables set above
if(alpine) {
  base_path <- paste0("/pl/active/jravilab/AGhosh/AMR_data/v2_data/", bug)
} else {
  base_path = "inst/misc"
}

mapping_file <- file.path(base_path, paste0("gene_", feature, "_mappings.tsv"))
count_file <- file.path(base_path, paste0(bug, "_gene_count.parquet"))
output_prefix <- file.path(base_path, paste0(bug, "_", feature))

# Load in the mapping data first
mapping_tbl <- read_tsv(mapping_file, col_types = cols()) %>%
  select(feature = target_name, gene = source) %>%
  distinct()

# And now load in the gene counts (gotta be counts, don't use binary please)
gene_counts <- read_parquet(count_file) %>%
  rename(gene = gene, genome = genome_id, count = value)

# Filter and map the genes to their COGs in the mapping table
mapped_counts <- gene_counts %>%
  inner_join(mapping_tbl, by = "gene") %>%
  group_by(genome, feature) %>%
  summarise(count = sum(count), .groups = "drop")

# Behold, a bog standard wide count matrix of your features of interest
count_matrix <- mapped_counts %>%
  pivot_wider(names_from = feature, values_from = count, values_fill = 0)

# Wow, a binary matrix too
binary_matrix <- count_matrix %>%
  mutate(across(-genome, ~ if_else(. > 0, 1L, 0L)))

# Save both matrices as compressed Parquet files for downstream matrix gen
write_parquet(count_matrix, paste0(output_prefix, "_counts.parquet"), compression = "zstd")
write_parquet(binary_matrix, paste0(output_prefix, "_binary.parquet"), compression = "zstd")
