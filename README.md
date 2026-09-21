
<!-- README.md is generated from README.Rmd. Please edit that file -->

# amRdata

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

**amRdata** is the first package in the [amR
suite](https://github.com/JRaviLab/amR) for antimicrobial resistance
(AMR) prediction. It takes a user‑provided species or taxon ID,
downloads the corresponding genomes and AST data from BV‑BRC, constructs
pangenomes, extracts features at multiple molecular scales, and prepares
a unified Parquet‑backed DuckDB file for downstream ML modeling in
**amRml**.

The workflow is comprised of 6 primary processes:

1.  BV‑BRC metadata (isolate metadata + AMR phenotypic labels) →
2.  BV-BRC genomes (sequence data) →
3.  Panaroo pangenome (genes, struct) →
4.  CD‑HIT protein clusters (proteins) →
5.  HMMER protein features (e.g., Pfam domains, COG membership, ARG
    homologs) →
6.  Database formatting

## Overview

**amRdata** includes functions to:

- Query and download bacterial genome data from BV-BRC
- Acquire paired antimicrobial susceptibility testing (AST) results
- Extract molecular features across scales:
  - Gene clusters (Panaroo pangenome analysis)
  - Protein clusters (CD-HIT sequence similarity)
  - Protein domains (Pfam annotations)
  - Structural variants (Panaroo pangenome rearrangements)
- Store all data in highly efficient Parquet and DuckDB formats

See the [package
vignette](https://jravilab.github.io/amRdata/articles/intro.html) for
detailed usage.

## Installation

``` r
# Install from GitHub
if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")

remotes::install_github("JRaviLab/amRdata")
```

## Quick start

``` r

library(amRdata)

# Step 1: Check BV-BRC data availability for bacteria of your choice
checkDataAvailability(
  user_bacs = c("Shigella flexneri", "Shigella sonnei", "Helicobacter pylori")
)

# Step 2: Download and prepare genomes with paired AST data from BV-BRC
prepareGenomes(
  user_bacs = c("Shigella flexneri"),
  base_dir  = "data/Shigella_flexneri",
  method    = "ftp",   # or "cli"
  verbose   = TRUE
)

# Step 3: Run full feature extraction (Panaroo → CD-HIT → HMMER → metadata cleaning)
runDataProcessing(
  duckdb_path = "data/Shigella_flexneri/Sfl.duckdb",
  output_path = "data/Shigella_flexneri",
  threads     = 16
)

# A final Parquet-backed DuckDB is created:
#   data/Shigella_flexneri/Sfl_parquet.duckdb

This contains data for feature presence/absence and counts across scales in 
genome by feature matrices, as well as all available sample metadata. These are
used in amRml, but can also be exported and used through:

# Optional — Step 5: Export your data in tabular format
### For basic stats and metadata
exportTables(
  duckdb_path = "data/Shigella_flexneri/Sfl.duckdb",
  export_formats = "csv"
)

### For all features and phenotypes after data processing
exportProcessedData(
  duckdb_path = "data/Shigella_flexneri/Sfl.duckdb",
  amr_phenotype_mode = "separate",
  export_formats = c("csv", "parquet", "xlsx"),
  export_sequences = TRUE
)
```

## Package features

### Data curation

1.  BV‑BRC data access **amRdata** calls the BV-BRC database through its
    API and FTP server by default, but can also utilize the BV‑BRC CLI
    (via Docker). It is an all-in-one method to fetch:

- Genome metadata
- AMR phenotype data
- Genome assemblies (`.fna`, `.faa`, `.gff`)

Users can query and fetch data through the functions:

    checkDataAvailability()
    prepareGenomes()

`prepareGenomes()` is a wrapper for the functions:

    retrieveMetadata()
    retrieveGenomes()

If the BV-BRC CLI is used instead of the default API, BV-BRC metadata is
cached automatically using `BiocFileCache` to speed up subsequent
CLI-based queries.

The package interfaces with BV-BRC (Bacterial and Viral Bioinformatics
Resource Center) to access bacterial genome sequences and antimicrobial
susceptibility testing data either using API/FTP or the BV-BRC CLI in a
Docker container for accessibility:

- Query isolate metadata with flexible filtering
- Download genome files (`.fna`, `.faa`, `.gff`)
- Retrieve AST results linking genotypes to phenotypes
- Apply quality control filters (assembly quality, metadata
  completeness)

### Feature extraction

Through a single user-friendly runner function `runDataProcessing()`,
features are extracted from downloaded genomes over multiple
complementary molecular scales:

#### 1. Gene clusters

Panaroo is executed inside a Docker container using `.runPanaroo()`.

Our pangenome creation approach:

- Allows end-to-end single pangenome runs
- Offers parallelized multi-batch pangenomes for large isolate sets
  (\>5,000 genomes)
  - Supports automated pangenome merging through `.mergePanaroo()`
- Generates gene presence/absence and count matrices per isolate
- Identifies structural variants (gene triplets indicating genome
  rearrangements)

#### 2. Protein clusters

CD-HIT is executed inside a Docker container using `.runCDHIT()`.

Our protein clustering approach:

- Clusters proteins across all isolates from BV-BRC `.faa` files
- Creates protein presence/absence and count matrices per isolate
- Saves cluster names and annotations

#### 3. HMMER features

HMMER is executed inside a Docker container using `.runHMMER`.

Our protein feature annotation approach:

- Automatically configures 4 HMMER databases for use
  - Pfam domain homology
  - Cluster of orthologous groups (COG) homology
  - Antimicrobial resistance gene (ARG) homology
  - Bacterial defense system homology
- Runs parallelized and containerized annotation
- Maps HMMER feature presence/absence and counts to genomes and proteins
- Provides multiple functional annotation layers

#### 4. Data cleaning and storage

Final data formatting and storage is executed using `cleanData()`.

Our final data storage script:

- Harmonizes drug names, classes, and countries in BV-BRC metadata
- Generates temporal bins to stratify analysis across time
- Summarizes AMR information across the dataset
- Writes all data into highly compressed data structures
  - **Parquet**: Binary, columnar storage for large matrices
    - These can be made human-readable by calling `arrow::read_parquet`
      or by using `exportTables()` and `exportProcessedData()`
  - **DuckDB**: SQL-queryable database for rapid filtering of linked
    Parquets

## Workflow example

An example of the minimal commands needed to download and process all
data and metadata for *Shigella flexneri* genomes with paired AST
metadata.

    library(amRdata)

    # 1. Download & filter genomes
    prepareGenomes(user_bacs = "Shigella flexneri")

    # 2. Run multi-scale feature extraction
    runDataProcessing(duckdb_path = "data/Shigella_flexneri/Sfl.duckdb")

    # This completes an amRdata analysis! 
    # Data are saved and ready for use or inspection.
    # For example, to view final metadata and features in human-readable format: 

    exportProcessedData(duckdb_path = "data/Shigella_flexneri/Sfl.duckdb", 
                        export_formats = c("tsv", "xlsx", "Parquet))

### A cautionary note about BV-BRC accession IDs

When reading exported tables, remember that BV-BRC accession IDs are
distinguished by trailing zeroes! This means 1280.10 and 1280.100 are
different genomes, but many programs will automatically truncate
trailing zeroes without warning a user.

Parquet is a highly efficient, quick file format, and it also avoids
this behavior! We recommend Parquet wherever possible, but as a
compressed, binary format, they need special handling to be
human-readable. See below.

### Reading binary Parquet files

    # To read Parquet files in R
    library(arrow)

    # To read recorded gene cluster counts recorded per isolate
    Sfl_gene_counts <- arrow::read_parquet("data/Shigella_flexneri/exported_data/gene_count.parquet")
      
    # To connect gene cluster IDs to their annotated names
    Sfl_gene_names <- arrow::read_parquet("data/Shigella_flexneri/exported_data/gene_names.parquet")

Some newer IDEs like Positron support loading human-readable Parquet
files by default.

## Data requirements

External dependencies (managed through Docker) <br>

- BV‑BRC CLI
- Panaroo
- CD‑HIT
- HMMER
- DuckDB
- Arrow (Parquet)

The user does not need to install these manually.

The package requires:

- An internet connection to access BV-BRC data and metadata
- A local Docker installation
  - Containers for internal tools are pulled automatically and do not
    require configuration
  - Make sure Docker is running before you start processing data!
- Sufficient storage for databases, downloaded files, and processed
  output (we recommend 50GB+)
  - These analyses can produce very large amounts of data
  - For taxa with \>5,000 genomes, disk space can easily exceed 100GB!
- Multicore processing and sufficient (16GB+) of RAM are highly
  recommended
  - Species with many isolates may run poorly or fail to complete on
    older hardware

### Output

Feature matrices dimensions depend on species:

- Rows: Number of isolates (typically \<10,000)
- Columns: Number of features (ballpark estimates)
  - Genes: 5,000-50,000
  - Proteins: 5,000-50,000
  - Pfam domains: 500-10,000
  - Structural variants: 1,000-10,000

### External dependencies

The package uses established bioinformatics tools:

- **Panaroo** (≥1.7.0): Pangenome analysis
- **CD-HIT** (≥4.8.1): Protein clustering
- **HMMER** (≥3.4): Protein feature annotation
- **Docker**: For BV-BRC CLI container

These are automatically managed through the Docker container.

## Performance

Processing times vary by species and isolate count:

- Data download: 0-1 hours

- Pangenome construction: 0-6 hours

- Protein clustering: 0-3 hours

- Protein feature annotation: 0-2 hours

- Total: 1-12 hours for a complete species analysis

- These numbers will all vary greatly based on isolate number, genome
  complexity, and available hardware.

- Parallelization significantly reduces processing time when multiple
  cores are available.

- If a `furrr`/`future::multisession` worker fails with
  `could not find function ".xxx"` for an internal amRdata helper, your
  installed copy of amRdata is stale relative to the source you’re
  editing. `future::multisession` workers are fresh R processes that
  resolve `amRdata` by loading the *installed* package from
  `.libPaths()` — they do not see changes made only via
  `devtools::load_all()` in your interactive session. Run
  `devtools::install()` (or
  `pkgbuild::compile_dll(); devtools::document(); devtools::install()`)
  before exercising any function that runs work via `future`/`furrr`, or
  temporarily set `future::plan(future::sequential)` while iterating
  with `load_all()` alone.

- `runDataProcessing()`/`runPanaroo2Duckdb()` default
  `panaroo_refind_mode` to `"off"` rather than Panaroo’s own default.
  Refinding recovers gene calls that annotation tools missed, but its
  search can take substantially longer (or in rare cases fail to
  complete within 6+ hours) when a genome carries a cluster of CDS with
  internal stop codons — a condition existing genome-quality metadata
  (CheckM completeness/contamination, consistency scores, quality flags)
  does not flag. This default trades some gene-recovery accuracy for
  predictable runtime; set `panaroo_refind_mode = "default"` to restore
  Panaroo’s normal behavior.

### Integration with amR suite

**amRdata** is designed to work seamlessly with other **amR** packages:

``` r
library(amRdata)
library(amRml)
library(amRviz)

# 1. Curate data
prepareGenomes("Shigella flexneri")
runDataProcessing("amRdata/data/Shigella_flexneri/Sfl.duckdb")

# 2. Train models
runMLmodels("amRdata/data/Shigella_flexneri/Sfl_parquet.duckdb")

# 3. Visualize ### To add
launch_dashboard()
```

## Related packages

- [amR](https://github.com/JRaviLab/amR): Suite metapackage
- [amRml](https://github.com/JRaviLab/amRml): ML for AMR prediction
- [amRviz](https://github.com/JRaviLab/amRviz): Interactive dashboard

## Citation

If you use `amRdata` in your research, please cite:

> Ghosh A^, Brenner EP^, Boyer EA, McKim AP, Vang CK, Wolfe EP, Mayer D,
> Lesiyon RL, Ravi J.
>
> amR: an R package suite to predict antimicrobial resistance in
> bacterial pathogens.
>
> bioRxiv. 2026. DOI:
> [10.64898/2026.07.10.734579](https://doi.org/10.64898/2026.07.10.734579).

^ Co-first authors

## Contributing

We welcome contributions! Please see [CONTRIBUTING.md](CONTRIBUTING.md)
for guidelines.

## Reporting issues

Report bugs and request features at:
<https://github.com/JRaviLab/amRml/issues>

## License

BSD 3-Clause License. See [LICENSE](LICENSE) for details.

## Contact

**Corresponding author**: Janani Ravi (<janani.ravi@cuanschutz.edu>)

**Lab website**: <https://jravilab.github.io>

## Code of conduct

Please note that **amRdata** is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
