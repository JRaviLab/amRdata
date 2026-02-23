computationalData <- function(duckdb_path, out_path){
    
con <- DBI::dbConnect(duckdb::duckdb(), duckdb_path)
on.exit(try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE), add = TRUE)
    
metadata <- DBI::dbReadTable(con, "metadata") |> 
    tibble::as_tibble() |>
    dplyr::mutate(
        `genome_drug.evidence` = dplyr::case_when(
          `genome_drug.laboratory_typing_method` %in%
            c("Disk diffusion", "MIC", "Broth dilution", "Agar dilution") ~ "Laboratory Method",
          `genome_drug.laboratory_typing_method` == "Computational Prediction" ~ "Computational Method",
          TRUE ~ `genome_drug.evidence`
        )
      )
    
 genomes_w_both_methods <- metadata |> 
    dplyr::select(genome_drug.genome_id, genome_drug.antibiotic, genome_drug.evidence) |> 
    dplyr::group_by(genome_drug.genome_id, genome_drug.antibiotic) |> 
    dplyr::count() |> 
    dplyr::filter(n == 2) |> 
    dplyr::ungroup() |> 
    dplyr::distinct(genome_drug.genome_id) 

DBI::dbWriteTable(con, "genomes_w_both_methods", genomes_w_both_methods, overwrite = TRUE) 
arrow::write_parquet(genomes_w_both_methods, file.path(out_path, "genomes_w_both_methods.parquet"), compression = "zstd", compression_level = 9, use_dictionary = TRUE)    
    
comp_only <- metadata |> 
    dplyr::select(genome_drug.genome_id, genome_drug.antibiotic, genome_drug.evidence) |> 
    dplyr::group_by(genome_drug.genome_id, genome_drug.antibiotic) |>
    dplyr::mutate(distinct_evidence = dplyr::n_distinct(genome_drug.evidence)) |>
    dplyr::ungroup() |>
    dplyr::filter(genome_drug.evidence == "Computational Method" & distinct_evidence == 1) |>
    dplyr::distinct(genome_drug.genome_id) 
    
DBI::dbWriteTable(con, "comp_only", comp_only, overwrite = TRUE) 
arrow::write_parquet(comp_only, file.path(out_path, "comp_only.parquet"), compression = "zstd", compression_level = 9, use_dictionary = TRUE)
   
}
