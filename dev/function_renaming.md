# amRdata function renaming proposal

## Rules

- Exported (external) names: lowerCamelCase, no leading dot.
- Internal names: leading dot plus lowerCamelCase. The dot is kept (it hides helpers from `ls()` and tab completion, and most internals already use it).
- Function names follow verbNoun (e.g. `.fetchBVBRCData`, `.buildPanarooGeneTable`, `.getHmmerVersion`). `is` predicates (e.g. `.isCompleteSet`) are accepted.
- Exceptions: the `.manifest*`, `.bfc*` and `.cli*` families keep their shared prefix instead of a leading verb, but use camelCase.
- Acronyms and initialisms are in uppercase, like the existing `.checkIDs`, `.extractAMRtable` and `.buildDBpath`: `ID`, `AMR`, `DB`, `QC`, `API`, `GFF`, `GTO`, `TSV`, `FASTA`, `ORB`, `HMMER`, `CDHIT`, `BVBRC`. `DuckDB` follows the product's own spelling. Words like `Metadata`, `Panaroo` stay in camelCase. We like `convertA2B`, but the verbNoun rule takes priority over `2`/`To` names: `.panaroo2duckdb` becomes `.writePanarooTables`, not `.panarooToDuckDB`.
- "—" means no change.

## All functions

| Existing | Suggested | Script | Type |
|---|---|---|---|
| `generateSummary` | `.generateSummary` | amRdataPlots.R | Internal |
| `generatePlots` | `.generatePlots` | amRdataPlots.R | Internal |
| `.writeMetadataReports` | `writeMetadataReports` | amRdataPlots.R | External |
| `.bvbrcApiReq` | `.sendBVBRCRequest` | bvbrc_api.R | Internal |
| `.bvbrcEnc` | `.encodeBVBRCQuery` | bvbrc_api.R | Internal |
| `.bvbrcApiFetch` | `.fetchBVBRCAPIData` | bvbrc_api.R | Internal |
| `.bvbrcChunk` | `.chunkBVBRCIDs` | bvbrc_api.R | Internal |
| `.bvbrcPrefixFill` | `.fillBVBRCPrefix` | bvbrc_api.R | Internal |
| `.extractAMRtableApi` | `.extractAMRtableFromAPI` | bvbrc_api.R | Internal |
| `.resolveGenomeIDsApi` | `.resolveGenomeIDsFromAPI` | bvbrc_api.R | Internal |
| `.extractGenomeDataApi` | `.extractGenomeDataFromAPI` | bvbrc_api.R | Internal |
| `.amr_lab_methods` | `.getLabMethods` | data_curation.R | Internal |
| `.create_amr_tagged_view` | `.createAMRTaggedView` | data_curation.R | Internal |
| `.ftps_download_one` | `.downloadFTPOne` | data_curation.R | Internal |
| `.ftps_download_two_pass` | `.downloadFTPTwoPass` | data_curation.R | Internal |
| `.fetchBVBRCdata` | `.fetchBVBRCData` | data_curation.R | Internal |
| `.parse_bvbrc_tsv` | `.parseBVBRCTSV` | data_curation.R | Internal |
| `.apply_metadata_qc` | `.applyMetadataQC` | data_curation.R | Internal |
| `.purge_genome_files` | `.purgeGenomeFiles` | data_curation.R | Internal |
| `.updateBVBRCdata` | — | data_curation.R | Internal |
| `.retrieveCustomQuery` | — | data_curation.R | Internal |
| `.generateDBname` | — | data_curation.R | Internal |
| `.amr_status` | `.reportStatus` | helpers.R | Internal |
| `.amr_progress_available` | `.isProgressAvailable` | helpers.R | Internal |
| `.amr_progress_handler` | `.getProgressHandler` | helpers.R | Internal |
| `.amr_with_progress` | `.evalWithProgress` | helpers.R | Internal |
| `.amr_download_file` | `.downloadFile` | helpers.R | Internal |
| `.amr_progressor` | `.createProgressor` | helpers.R | Internal |
| `.amr_progress_step` | `.runProgressStep` | helpers.R | Internal |
| `.amr_dataset_paths` | `.getDatasetPaths` | helpers.R | Internal |
| `.amr_paths_from_duckdb` | `.getDuckDBPaths` | helpers.R | Internal |
| `.amr_paths_from_dataset_db` | `.getDatasetDBPaths` | helpers.R | Internal |
| `.amr_connect_dataset_db` | `.connectDatasetDB` | helpers.R | Internal |
| `.buildDBpath` | — | helpers.R | Internal |
| `.retrieveQueryIDs` | — | data_curation.R | Internal |
| `.extractAMRtable` | — | data_curation.R | Internal |
| `.extractGenomeData` | — | data_curation.R | Internal |
| `retrieveMetadata` | — | data_curation.R | External |
| `.filterGenomes` | — | data_curation.R | Internal |
| `.is_complete_set` | `.isCompleteSet` | data_curation.R | Internal |
| `.list_complete` | `.listCompleteGenomes` | data_curation.R | Internal |
| `.cli_dump_fastas_gto_chunk` | `.cliDumpFASTAsGTOChunk` | data_curation.R | Internal |
| `.cli_export_gff_chunk` | `.cliExportGFFChunk` | data_curation.R | Internal |
| `retrieveGenomes` | — | data_curation.R | External |
| `genomeList` | `.buildGenomeList` | data_curation.R | Internal |
| `prepareGenomes` | — | data_curation.R | External |
| `exportTables` | — | data_curation.R | External |
| `checkDataAvailability` | — | data_curation.R | External |
| `removeLocalDatabases` | — | data_processing.R | External |
| `.processPanaroo` | `.launchPanaroo` | data_processing.R | Internal |
| `.runPanaroo` | `.orchestratePanaroo` | data_processing.R | Internal |
| `.mergePanaroo` | — | data_processing.R | Internal |
| `.panaroo2geneTable` | `.buildPanarooGeneTable` | data_processing.R | Internal |
| `.panaroo2geneNames` | `.buildPanarooGeneNames` | data_processing.R | Internal |
| `.panaroo2StructTable` | `.buildPanarooStructTable` | data_processing.R | Internal |
| `.panaroo2OtherTables` | `.buildPanarooOtherTable` | data_processing.R | Internal |
| `.panaroo2duckdb` | `.writePanarooTables` | data_processing.R | Internal |
| `.runCDHIT` | `.launchCDHIT` | data_processing.R | Internal |
| `runPanaroo2Duckdb` | `extractPangenomeFeatures` | data_processing.R | External |
| `.parseProteinClusters` | — | data_processing.R | Internal |
| `.extractMembersInClusters` | — | data_processing.R | Internal |
| `.buildProtMatrices` | `.buildProteinMatrix` | data_processing.R | Internal |
| `buildMatrices` | remove (back-compat wrapper for `.buildProteinMatrix`) | data_processing.R | Internal |
| `.clusterNames` | `.getClusterName` | data_processing.R | Internal |
| `CDHIT2duckdb` | `extractProteinFeatures` | data_processing.R | External |
| `.prepareHmmerDatabases` | — | data_processing.R | Internal |
| `.runHmmerJob` | `.launchHMMERJob` | data_processing.R | Internal |
| `.runHMMER` | `.searchHMMERDBs` | data_processing.R | Internal |
| `.proteinAnnotations2Duckdb` | `.writeProteinAnnotations` | data_processing.R | Internal |
| `.defenseHMMER` | `.searchDefenseCas` | data_processing.R | Internal |
| `HMMER2duckdb` | `annotateProteins` | data_processing.R | External |
| `cleanMetaData` | `cleanMetadata` | data_processing.R | External |
| `cleanData` | `cleanGenomeData` | data_processing.R | External |
| `runDataProcessing` | `processGenomeFeatures` | data_processing.R | External |
| `exportProcessedData` | — | data_processing.R | External |
| `removeLocalFiles` | — | data_processing.R | External |
| `buildDyadFeatureMap` | — | feature_to_head.R | External |
| `.resolve_workers` | `.resolveWorkers` | helpers.R | Internal |
| `.amr_set_future_plan` | `.setFuturePlan` | helpers.R | Internal |
| `.bvbrcFutureMap` | — | helpers.R | Internal |
| `.id_checker` | `.checkIDs` | helpers.R | Internal |
| `.preserve_export_id_text` | `.preserveExportIDText` | helpers.R | Internal |
| `.amr_bfc` | `.bfc` | helpers.R | Internal |
| `.amr_bfc_find` | `.bfcFind` | helpers.R | Internal |
| `.amr_bfc_bvbrc_path` | `.bfcBVBRCPath` | helpers.R | Internal |
| `.amr_bfc_remove_bvbrc` | `.bfcRemoveBVBRC` | helpers.R | Internal |
| `.amr_bfc_hmmer_rname` | `.bfcHMMERRname` | helpers.R | Internal |
| `.amr_bfc_bvbrc_rname` | `.bfcBVBRCRname` | helpers.R | Internal |
| `.amr_bfc_register_hmmer` | `.bfcRegisterHMMER` | helpers.R | Internal |
| `.amr_bfc_hmmer_resources` | `.bfcHMMERResources` | helpers.R | Internal |
| `.amr_bfc_register_local` | `.bfcRegisterLocal` | helpers.R | Internal |
| `.amr_processing_datasets` | `.listProcessingDatasets` | helpers.R | Internal |
| `.amr_cleanup_datasets` | `.listCleanupDatasets` | helpers.R | Internal |
| `.amr_select_processing_dataset` | `.selectProcessingDataset` | helpers.R | Internal |
| `.docker_path` | `.getDockerPath` | helpers.R | Internal |
| `.pick_shell` | `.pickShell` | helpers.R | Internal |
| `.strip_fasta_preamble` | `.stripFASTAPreamble` | helpers.R | Internal |
| `.sanitize_gff` | `.sanitizeGFF` | helpers.R | Internal |
| `.checkDataPerTaxon` | — | helpers.R | Internal |
| `.manifest_file_info` | `.manifestFileInfo` | helpers.R | Internal |
| `.manifest_git_info` | `.manifestGitInfo` | helpers.R | Internal |
| `.manifest_package_versions` | `.manifestPackageVersions` | helpers.R | Internal |
| `.manifest_run_id` | `.manifestRunID` | helpers.R | Internal |
| `.manifest_start` | `.manifestStart` | helpers.R | Internal |
| `.amr_bfc_manifest_rname` | `.bfcManifestRname` | helpers.R | Internal |
| `.manifest_artifact` | `.manifestArtifact` | helpers.R | Internal |
| `.manifest_validate` | `.manifestValidate` | helpers.R | Internal |
| `.manifest_stage` | `.manifestStage` | helpers.R | Internal |
| `.manifest_event` | `.manifestEvent` | helpers.R | Internal |
| `.manifest_finish` | `.manifestFinish` | helpers.R | Internal |
| `.log_write` | `.writeLog` | helpers.R | Internal |
| `.log_tool_output` | `.logTools` | helpers.R | Internal |
| `.log_or_message` | `.logOrMessage` | helpers.R | Internal |
| `.manifest_prior_stage` | `.manifestPriorStage` | helpers.R | Internal |
| `.resume_plan` | `.resumePlan` | helpers.R | Internal |
| `.manifest_find_latest` | `.manifestFindLatest` | helpers.R | Internal |
| `.manifest_resume` | `.manifestResume` | helpers.R | Internal |
| `.to_container` | `.convert2ContainerPath` | helpers.R | Internal |
| `.stripPseudogeneGFFs` | `.stripPseudogeneGFF` | helpers.R | Internal |
| `.exportDyadAnnotations` | — | helpers.R | Internal |
| `.amr_resolve_export_orb` | `.resolveExportORB` | helpers.R | Internal |
| `.readOrbMetadata` | `.readORBMetadata` | helpers.R | Internal |
| `.isValidHmmFile` | `.isValidHMMFile` | helpers.R | Internal |
| `.parse_hmmer_profiles` | `.parseHMMERProfiles` | helpers.R | Internal |
| `.parseHMMEROutput` | — | helpers.R | Internal |
| `.write_compressed_parquet` | `.writeCompressedParquet` | helpers.R | Internal |
| `.defaultHmmerDbDir` | `.getHMMERDBDir` | helpers.R | Internal |
| `.hmmer_version` | `.getHMMERVersion` | helpers.R | Internal |
| `meta_palette` | `.getMetaPalette` | utils_colors.R | Internal |


