library(dplyr)
library(readr)
library(stringr)
library(arrow)
library(fs)
library(purrr)

# Set directories
top_feat_dir <- "results/tsv"
species <- "Aba"  # Change as needed
seq_dir <- file.path("/pl/active/jravilab/AGhosh/AMR_data/v2_data", species)
output_dir <- "hmmer_inputs"
#dir_create(output_dir)

# Load all top feature files
top_files <- dir_ls(top_feat_dir, regexp = "_top_features\\.tsv$")

# Step 1: Aggregate all unique feature IDs
all_features <- top_files %>%
  map_dfr(~ read_tsv(.x, show_col_types = FALSE) %>%
            select(Variable) %>%
            mutate(source_file = path_file(.x))) %>%
  distinct()

# Step 2: Load gene and protein sequences
gene_seqs <- read_parquet(file.path(seq_dir, paste0(species, "_gene_seqs.parquet")))
protein_seqs <- read_parquet(file.path(seq_dir, paste0(species, "_protein_seqs.parquet")))

# Join features to sequences
gene_matches <- all_features %>%
  inner_join(gene_seqs, by = c("Variable" = "name")) %>%
  mutate(type = "gene")

protein_matches <- all_features %>%
  inner_join(protein_seqs, by = c("Variable" = "name")) %>%
  mutate(type = "protein")

# Combine and save mapping table
mapping_table <- bind_rows(gene_matches, protein_matches)
write_tsv(mapping_table, file.path(output_dir, paste0(species, "_feature_sequence_mapping.tsv")))

# Write FASTA files for HMMer
write_lines(
  paste0("> ", gene_matches$Variable, "\n", gene_matches$sequence),
  file.path(output_dir, paste0(species, "_genes_for_hmmer.fasta"))
)

write_lines(
  paste0("> ", protein_matches$Variable, "\n", protein_matches$sequence),
  file.path(output_dir, paste0(species, "_proteins_for_hmmer.fasta"))
)

message("Mapping table and FASTA files created for HMMer input.")


### New try!
library(dplyr)
library(readr)
library(arrow)
library(fs)

# --- Config ---
species <- "Aba"
seq_dir <- file.path("/pl/active/jravilab/AGhosh/AMR_data/v2_data", species)
output_dir <- "hmmer_inputs"
dir_create(output_dir)

# --- Load all gene and protein sequences ---
gene_seqs <- read_parquet(file.path(seq_dir, paste0(species, "_gene_seqs.parquet"))) %>%
  mutate(type = "gene")

protein_seqs <- read_parquet(file.path(seq_dir, paste0(species, "_protein_seqs.parquet"))) %>%
  mutate(type = "protein")

# --- Combine and save mapping table ---
mapping_table <- bind_rows(gene_seqs, protein_seqs) %>%
  rename(Variable = name)

write_tsv(mapping_table, file.path(output_dir, paste0(species, "_all_feature_sequence_mapping.tsv")))

# --- Write FASTA files for HMMer ---
write_lines(
  paste0("> ", gene_seqs$name, "\n", gene_seqs$sequence),
  file.path(output_dir, paste0(species, "_genes_for_hmmer.fasta"))
)

write_lines(
  paste0("> ", protein_seqs$name, "\n", protein_seqs$sequence),
  file.path(output_dir, paste0(species, "_proteins_for_hmmer.fasta"))
)

message("All sequences exported for HMMer input.")

### Advanced try
# Generate job lists for ARG and COG models across species
library(dplyr)
library(readr)
library(stringr)
library(arrow)
library(fs)
library(purrr)

species_list <- c("Efa", "Sau", "Kpn", "Aba", "Pae", "Esp.")
chunk_count <- 32
input_dir <- "hmmer_inputs"
db_map <- c(
  COG = "COG_database2024.hmm",
  ARG = "NCBI_ResFinder.hmm"
)

dir_create(input_dir)

for (species in species_list) {
  seq_dir <- file.path("/pl/active/jravilab/AGhosh/AMR_data/v2_data", species)

  gene_seqs <- read_parquet(file.path(seq_dir, paste0(species, "_gene_seqs.parquet"))) %>%
    mutate(type = "gene")
  protein_seqs <- read_parquet(file.path(seq_dir, paste0(species, "_protein_seqs.parquet"))) %>%
    mutate(type = "protein")

  # Write full FASTA files
  write_lines(paste0("> ", gene_seqs$name, "\n", gene_seqs$sequence),
              file.path(input_dir, paste0(species, "_genes_for_hmmer.fasta")))
  write_lines(paste0("> ", protein_seqs$name, "\n", protein_seqs$sequence),
              file.path(input_dir, paste0(species, "_proteins_for_hmmer.fasta")))

  # Chunking function
  split_fasta <- function(seqs, prefix) {
    records <- paste0("> ", seqs$name, "\n", seqs$sequence)
    chunk_size <- ceiling(length(records) / chunk_count)
    chunks <- split(records, ceiling(seq_along(records) / chunk_size))

    walk2(chunks, seq_along(chunks), function(chunk, i) {
      chunk_path <- file.path(input_dir, sprintf("%s_chunk_%02d.fasta", prefix, i))
      write_lines(chunk, chunk_path)
    })
  }

  split_fasta(gene_seqs, paste0(species, "_genes"))
  split_fasta(protein_seqs, paste0(species, "_proteins"))

  # Generate job list for ARG and COG
  job_list <- expand.grid(
    type = c("genes", "proteins"),
    chunk = sprintf("%02d", 1:chunk_count),
    db = names(db_map),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      JOB_NAME = paste0(type, "_chunk_", chunk, "_", db),
      FASTA = paste0(species, "_", type, "_chunk_", chunk, ".fasta"),
      DB = db_map[db]
    ) %>%
    select(JOB_NAME, FASTA, DB)

  write_tsv(job_list, file.path(input_dir, paste0(species, "_hmmer_jobs_COG_ARG.txt")))
}
