# amRdata function renaming proposal

## Rules

- Exported (external) names: lowerCamelCase, no leading dot.
- Internal names: leading dot plus lowerCamelCase. The dot is kept (it hides helpers from `ls()` and tab completion, and most internals already use it).
- Function names follow verbNoun (e.g. `.fetchBvbrcData`, `.buildPanarooGeneTable`, `.getHmmerVersion`). `is` predicates (e.g. `.isCompleteSet`) are accepted.
- Exceptions: the `.manifest*`, `.bfc*` and `.cli*` families keep their shared prefix instead of a leading verb, but use camelCase.
- Acronyms: `Hmmer`, `Cdhit`, `Duckdb`, `Metadata`, and `ToX` instead of `2x`.
- "—" means no change.

## All functions

| Existing | Suggested | Script | Type |
|---|---|---|---|
| `generateSummary` | — | amRdataPlots.R | External |
| `generatePlots` | — | amRdataPlots.R | External |
| `.bvbrcApiReq` | `.sendBvbrcRequest` | bvbrc_api.R | Internal |
| `.bvbrcEnc` | `.encodeBvbrcQuery` | bvbrc_api.R | Internal |
| `.bvbrcApiFetch` | `.fetchBvbrcData` | bvbrc_api.R | Internal |
| `.bvbrcChunk` | `.chunkBvbrcIDs` | bvbrc_api.R | Internal |
| `.bvbrcPrefixFill` | `.fillBvbrcPrefix` | bvbrc_api.R | Internal |
| `.extractAMRtableApi` | `.extractAMRtableFromApi` | bvbrc_api.R | Internal |
| `.resolveGenomeIDsApi` | `.resolveGenomeIDsFromApi` | bvbrc_api.R | Internal |
| `.extractGenomeDataApi` | `.extractGenomeDataFromApi` | bvbrc_api.R | Internal |
| `.amr_lab_methods` | `.getLabMethods` | data_curation.R | Internal |
| `.create_amr_tagged_view` | `.createAMRTaggedView` | data_curation.R | Internal |
| `.ftps_download_one` | `.downloadFtpsOne` | data_curation.R | Internal |
| `.ftps_download_two_pass` | `.downloadFtpsTwoPass` | data_curation.R | Internal |
| `.fetchBVBRCdata` | — | data_curation.R | Internal |
| `.parse_bvbrc_tsv` | `.parseBvbrcTsv` | data_curation.R | Internal |
| `.apply_metadata_qc` | `.applyMetadataQC` | data_curation.R | Internal |
| `.purge_genome_files` | `.purgeGenomeFiles` | data_curation.R | Internal |
| `.ensure_bvbrc_cache` | `.ensureBvbrcCache` | data_curation.R | Internal |
| `.updateBVBRCdata` | - | data_curation.R | Internal |
| `.retrieveCustomQuery` | — | data_curation.R | Internal |
| `.generateDBname` | - | data_curation.R | Internal |
| `.buildDBpath` | - | data_curation.R | Internal |
| `.retrieveQueryIDs` | — | data_curation.R | Internal |
| `.extractAMRtable` | - | data_curation.R | Internal |
| `.extractGenomeData` | — | data_curation.R | Internal |
| `retrieveMetadata` | — | data_curation.R | External |
| `.filterGenomes` | — | data_curation.R | Internal |
| `.is_complete_set` | `.isCompleteSet` | data_curation.R | Internal |
| `.list_complete` | `.listCompleteGenomes` | data_curation.R | Internal |
| `.cli_dump_fastas_gto_chunk` | `.cliDumpFastasGtoChunk` | data_curation.R | Internal |
| `.cli_export_gff_chunk` | `.cliExportGffChunk` | data_curation.R | Internal |
| `retrieveGenomes` | — | data_curation.R | External |
| `genomeList` | `.buildGenomeList` | data_curation.R | Internal |
| `prepareGenomes` | — | data_curation.R | External |
| `exportTables` | - | data_curation.R | External |
| `checkDataAvailability` | — | data_curation.R | External |
| `clearHMMERdatabases` | `clearHmmerDatabases` | data_processing.R | External |
| `.processPanaroo` | `.launchPanaroo` | data_processing.R | Internal |
| `.runPanaroo` | `.orchestratePanaroo` | data_processing.R | Internal |
| `.mergePanaroo` | - | data_processing.R | Internal |
| `.panaroo2geneTable` | `.buildPanarooGeneTable` | data_processing.R | Internal |
| `.panaroo2geneNames` | `.buildPanarooGeneNames` | data_processing.R | Internal |
| `.panaroo2StructTable` | `.buildPanarooStructTable` | data_processing.R | Internal |
| `.panaroo2OtherTables` | `.buildPanarooOtherTables` | data_processing.R | Internal |
| `.panaroo2duckdb` | `.writePanarooTables` | data_processing.R | Internal |
| `.runCDHIT` | `.launchCdhit` | data_processing.R | Internal |
| `runPanaroo2Duckdb` | `extractPangenomeFeatures` | data_processing.R | External |
| `.parseProteinClusters` | — | data_processing.R | Internal |
| `.extractMembersInClusters` | — | data_processing.R | Internal |
| `.buildProtMatrices` | `.buildProteinMatrices` | data_processing.R | Internal |
| `.clusterNames` | `.getClusterNames` | data_processing.R | Internal |
| `CDHIT2duckdb` | `extractProteinFeatures` | data_processing.R | External |
| `.prepareHmmerDatabases` | — | data_processing.R | Internal |
| `.runHmmerJob` | `.launchHmmerJob` | data_processing.R | Internal |
| `.runHMMER` | `.searchHmmerDatabases` | data_processing.R | Internal |
| `.proteinAnnotations2Duckdb` | `.writeProteinAnnotations` | data_processing.R | Internal |
| `.defenseHMMER` | `.searchDefenseCas` | data_processing.R | Internal |
| `cleanMetaData` | `cleanMetadata` | data_processing.R | External |
| `cleanData` | `cleanGenomeData` | data_processing.R | External |
| `runDataProcessing` | `processGenomeFeatures` | data_processing.R | External |
| `exportProcessedData` | — | data_processing.R | External |
| `buildDyadFeatureMap` | — | feature_to_head.R | External |
| `.resolve_workers` | `.resolveWorkers` | helpers.R | Internal |
| `.amr_set_future_plan` | `.setFuturePlan` | helpers.R | Internal |
| `.bvbrcFutureMap` | — | helpers.R | Internal |
| `.id_checker` | `.checkIDs` | helpers.R | Internal |
| `.preserve_export_id_text` | `.preserveExportIdText` | helpers.R | Internal |
| `.amr_bfc` | `.bfc` | helpers.R | Internal |
| `.amr_bfc_find` | `.bfcFind` | helpers.R | Internal |
| `.amr_bfc_bvbrc_path` | `.bfcBvbrcPath` | helpers.R | Internal |
| `.amr_bfc_hmmer_rname` | `.bfcHmmerRname` | helpers.R | Internal |
| `.amr_bfc_register_hmmer` | `.bfcRegisterHmmer` | helpers.R | Internal |
| `.amr_bfc_hmmer_resources` | `.bfcHmmerResources` | helpers.R | Internal |
| `.amr_bfc_register_local` | `.bfcRegisterLocal` | helpers.R | Internal |
| `.docker_path` | `.getDockerPath` | helpers.R | Internal |
| `.pick_shell` | `.pickShell` | helpers.R | Internal |
| `.strip_fasta_preamble` | `.stripFastaPreamble` | helpers.R | Internal |
| `.sanitize_gff` | `.sanitizeGff` | helpers.R | Internal |
| `.checkDataPerTaxon` | — | helpers.R | Internal |
| `.manifest_file_info` | `.manifestFileInfo` | helpers.R | Internal |
| `.manifest_git_info` | `.manifestGitInfo` | helpers.R | Internal |
| `.manifest_package_versions` | `.manifestPackageVersions` | helpers.R | Internal |
| `.manifest_run_id` | `.manifestRunId` | helpers.R | Internal |
| `.manifest_start` | `.manifestStart` | helpers.R | Internal |
| `.amr_bfc_manifest_rname` | `.bfcManifestRname` | helpers.R | Internal |
| `.manifest_artifact` | `.manifestArtifact` | helpers.R | Internal |
| `.manifest_migrate_legacy` | `.manifestMigrateLegacy` | helpers.R | Internal |
| `.manifest_validate` | `.manifestValidate` | helpers.R | Internal |
| `.manifest_stage` | `.manifestStage` | helpers.R | Internal |
| `.manifest_event` | `.manifestEvent` | helpers.R | Internal |
| `.manifest_finish` | `.manifestFinish` | helpers.R | Internal |
| `.log_write` | `.writeLog` | helpers.R | Internal |
| `.manifest_prior_stage` | `.manifestPriorStage` | helpers.R | Internal |
| `.resume_plan` | `.resumePlan` | helpers.R | Internal |
| `.manifest_find_latest` | `.manifestFindLatest` | helpers.R | Internal |
| `.manifest_resume` | `.manifestResume` | helpers.R | Internal |
| `.to_container` | `.convertToContainerPath` | helpers.R | Internal |
| `.stripPseudogeneGFFs` | `.stripPseudogeneGff` | helpers.R | Internal |
| `.exportDyadAnnotations` | — | helpers.R | Internal |
| `.isValidHmmFile` | — | helpers.R | Internal |
| `.parse_hmmer_profiles` | `.parseHmmerProfiles` | helpers.R | Internal |
| `.parseHMMEROutput` | `.parseHmmerOutput` | helpers.R | Internal |
| `.write_compressed_parquet` | `.writeCompressedParquet` | helpers.R | Internal |
| `.defaultHmmerDbDir` | `.getHmmerDbDir` | helpers.R | Internal |
| `.hmmer_version` | `.getHmmerVersion` | helpers.R | Internal |
| `meta_palette` | `.getMetaPalette` | utils_colors.R | Internal |


