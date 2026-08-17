# BV-BRC Data API — Feasibility & Limitations for amRdata

**Status:** evaluation notes (not yet implemented)
**Date:** 2026-07-24
**Author:** Emily Boyer
**Scope:** Can we replace the current Docker/`p3-*` CLI download path with direct
calls to the BV-BRC Data API, to fix the download bottleneck and let users pull
genome metadata + AMR phenotype data together, choosing their own columns?

All numbers below were measured with live `curl` calls against
`https://www.bv-brc.org/api/` on the date above. No authentication is required
for public data.

---

## 1. Bottom line

**The API is viable and removes the bottleneck. We should call it directly from R
(`httr2`), and NOT adopt the `bvbrc` Python package.**

- The API is a plain REST/Solr service. Everything we need — column selection,
  counting without downloading, faceted subsampling — is a URL query.
- The Python package (`bvbrc` v0.2.1, Beta, 1 maintainer, 9 commits) is a thin
  wrapper over this same API. It offers no capability R can't reach with `httr2`,
  and it does **not** solve the one hard part (the 25,000-row cap — see §4). It
  would add a `reticulate` + Python 3.9–3.11 dependency for no net gain. Use its
  client design as a reference, not as a dependency.

---

## 2. What our use case needs → feasibility

Derived from the AMR team meeting notes.

| Requirement | Feasible | Mechanism |
|---|---|---|
| User specifies **which columns** to return | ✅ | `select(field1,field2,...)` |
| Count matches **without downloading** ("-n" option) | ✅ | `Content-Range` response header — instant |
| Pick **species/pathogen from a list** (not typed) | ✅ | Solr faceting on `species` (genome collection) |
| Subsample "**500 genomes with most coverage** for drugs of interest" | ✅ | Facet `genome_id` filtered by `in(antibiotic,(...))`, ranked by count |
| Pull **QC-passing snapshot** (Good genomes), relevant columns → Zenodo | ✅ | `eq(genome_quality,Good)` + `select(...)` + keyset paging (§4) |
| Get metadata + AMR phenotype **together** | ⚠️ partial | Two collections (`genome`, `genome_amr`); no server-side join. But each call is sub-second, so it's 2 fast pulls + a local join instead of 2 slow Docker pipelines |

The "together" goal isn't literally a single call, but it stops being a
bottleneck: today's pain is Docker container startup + serial TSV parsing, not
the join.

---

## 3. Verified capabilities & scale

| Metric | Value |
|---|---|
| Total `genome_amr` rows (whole DB) | 17,266,649 |
| Total `genome` records | 16,885,561 |
| "Good" quality genomes | 8,117,447 |
| Throughput | 25,000 rows / ~1.9 s / ~2.2 MB (4 columns) |
| Column selection (`select`) | works |
| Count-only (`Content-Range`) | works, instant |
| Faceting (species list; genome ranking by drug coverage) | works |
| Auth for public data | none required |

**Query interface:** RQL over HTTP GET/POST, e.g.
`.../genome_amr/?and(eq(genome_name,Staphylococcus aureus),in(antibiotic,(ciprofloxacin,gentamicin)))&select(genome_id,antibiotic,resistant_phenotype)&sort(+id)&limit(25000)`
(URL-encode field names and values individually). POST the query body for large
`in(...)` lists to avoid URL-length limits.

---

## 4. The key limitation: the 25,000-row cap — and how to get around it

**Two hard ceilings, both = 25,000:**

1. **Per-request cap.** `limit(30000)` returns only 25,000 rows.
2. **Offset ceiling.** `limit(count, start)` with `start >= 25000` returns
   **HTTP 400**. So you cannot walk a large result set by increasing the offset.

`cursorMark` (Solr's deep-scroll) is **not reachable through the RQL layer** —
passing `&cursorMark=*` errors with `undefined field: "cursorMark"`.

### Answer to "we don't know how much we'll need to download at a given time"

Use a **two-step pattern** that scales to any size without knowing it in advance:

**Step 1 — Count first (free, instant).** Read the `Content-Range` header on a
`limit(1)` request. `items 0-1/590284` means 590,284 matching rows. Now you know
the exact size before downloading anything. This directly resolves the "we don't
know how much" problem — you always find out up front, for zero cost.

**Step 2 — If the count > 25,000, use keyset (seek) pagination.** Instead of an
offset, sort by the unique `id` and, each page, ask for rows *after* the last id
you saw. `start` is always 0, so the offset ceiling never applies. Loop until a
page returns fewer than 25,000 rows. This handles arbitrary/unknown totals.

```
# Pseudocode — retrieve ALL rows for a query, any size
total   <- GET .../genome_amr/?<query>&limit(1)   # read Content-Range -> N
rows    <- []
last_id <- ""                                     # empty sorts before all ids
repeat:
    page <- GET .../genome_amr/?and(<query>, gt(id, last_id))
                 &select(<cols>)&sort(+id)&limit(25000)
    if last_id != "":            # gt(id,...) is INCLUSIVE on this API
        drop page[0] if page[0].id == last_id   # dedupe boundary row
    rows    <- rows + page
    last_id <- page[last].id
    until length(page) < 25000
```

**Gotcha (verified):** `gt(id, X)` behaves as `>=` here — the boundary row
reappears as the first row of the next page. Dedupe on `id` when stitching pages,
or drop the first row of every page after the first.

**Cost example:** the full 17.3M-row `genome_amr` collection is ~692 pages of
25k. At ~1.9 s/page that's ~22 min sequential for the entire AMR table — or a few
minutes with parallel id-partitioning (see below). And we never need the whole
thing; a QC-filtered, column-selected subset is far smaller. Contrast with the
current per-500-genome Docker-container approach.

### Alternative: partition the query
If a single query is huge, split it by a natural key (species, taxon, or year)
so each sub-query is < 25k, and page each partition. Keyset paging is simpler and
preferred; partitioning is a fallback when one key range dominates.

### Looping vs. parallelization (both measured)

**Looping is the keyset walk itself, and it is inherently sequential *within* a
result set:** page N needs the last `id` from page N-1, so a single keyset walk
cannot be parallelized.

**Parallelization works by partitioning the `id` space** into disjoint ranges and
giving each worker its own range. It does **not** work by offset — the 25k offset
ceiling (§4) rules out `limit(25000, 25000)`-style parallel offsets entirely.

Because `id` is a UUID (uniformly distributed hex), splitting on the first hex
character gives 16 disjoint, naturally **balanced** buckets — no knowledge of the
data distribution required. Verified bucket sizes:

- rows with `id` in `[0,1)`: 1,078,515
- rows with `id` in `[1,2)`: 1,078,740
- × 16 buckets ≈ 17.3M total ✓

Each bucket is bounded with `and(ge(id,<lo>),lt(id,<hi>))` and is still >25k, so a
worker **keyset-paginates within its bucket**. Two-level structure:

```
partition id space -> [0,1), [1,2), ... [f,g)      # 16 disjoint buckets
  |-- worker pool (BiocParallel / future) --|
      each worker: keyset-paginate its bucket(s) with gt(id,last) + sort(+id)
```

**Measured speedup:** 4 workers pulling 25k-row pages from 4 buckets took **5.1 s
vs 18.5 s sequential (~3.6×), with no throttling observed.** This maps directly
onto the parallel backend the package already uses (`BiocParallel` / `future`).

**Be a good citizen.** Rate limits are undocumented, so keep concurrency modest
(≈4–8 workers, not 16+), reuse connections, and add retry-with-backoff on
transient failures (HTTP 5xx) — the same robustness gap that bites the current
Docker path. Finer partitioning (2 hex chars = 256 buckets) is available if you
want a work queue that load-balances across a fixed worker pool.

---

## 5. Other limitations / open items

- **`genome` and `genome_amr` are separate collections.** No server-side join;
  join locally on `genome_id`. Fine, since each pull is fast.
- **Distinct-genome counts** (e.g. "how many genomes actually have AMR data" — the
  true Zenodo snapshot size) need Solr's native `json.facet`/`unique()`, which the
  RQL facet layer did not return cleanly. Obtainable via the native Solr interface
  (same host, `Accept: application/solr+json`); falls out of the same work as deep
  scrolling. Not yet measured. `unique()` is an HLL estimate, not exact.
- **Genome sequence files** (`.fna`/`.faa`/`.gff`) are a separate concern — the
  Data API returns metadata/phenotype records, not assembly files. Those still
  come from FTP/CLI. This evaluation covers metadata + AMR phenotype only.
- **Rate limits** are not documented; be a good citizen (reuse connections, avoid
  hammering, prefer count-first over speculative pulls).

---

## 6. Recommendation

1. Build a small R module (`httr2`) with:
   - a query builder (RQL: `eq`/`in`/`gt`/`select`/`sort`/`limit`),
   - a `count()` helper reading `Content-Range`,
   - a `fetch_all()` that keyset-paginates past 25k with boundary dedupe, with
     optional id-space partitioning across a parallel worker pool for bulk pulls,
   - faceting helpers for the "pick a species" and "rank genomes by drug
     coverage" subsampling steps.
2. Benchmark it against the current Docker `.extractAMRtable()` /
   `.extractGenomeData()` path.
3. Do **not** take on the Python package as a dependency.

---

## 7. Prototype & benchmark

A runnable prototype lives at [`dev/bvbrc_api_prototype.R`](../dev/bvbrc_api_prototype.R)
(standalone; not wired into the package). It implements `bvbrc_count()`,
`bvbrc_fetch_all()` (keyset + optional parallel id-partitioning), and
`bvbrc_rank_genomes_by_drug_coverage()`. Run `bvbrc_demo()` after sourcing.

**Measured (all S. aureus AMR rows):**

| Operation | Result |
|---|---|
| `bvbrc_count()` | 590,284 rows — instant (header only) |
| `bvbrc_fetch_all()` sequential | 590,284 rows in **72 s** |
| `bvbrc_fetch_all(parallel=TRUE, workers=6)` | 590,284 rows in **27 s (~2.6×)**, identical row count |
| keyset correctness (small set, tiny page size) | rows == count, all ids unique, no boundary duplication |

Parallel returns the exact same row count as sequential — the id-space
partitioning is complete and non-overlapping.

---

## 8. Other amRdata use cases the API covers

The API is not just for AMR phenotype — it reaches most of what the package
currently shells out to Docker (`p3-*`) for. All verified live:

| Need (current source) | API collection | Verified |
|---|---|---|
| AMR phenotype (`genome_amr` via `p3-get-genome-drugs`) | `genome_amr` | ✅ |
| Genome metadata (via `p3-get-genome-data`) | `genome` | ✅ |
| Gene/annotation content — `.PATRIC.gff` (via `p3-dump-genomes`) | `genome_feature` | ✅ 2,583 CDS features w/ `patric_id`, `product`, coords, strand, `aa_sequence_md5` |
| **Protein FASTA — `.PATRIC.faa`** (via `p3-dump-genomes`) | `genome_feature` → `feature_sequence` | ✅ md5 → actual AA sequence returned |
| Contig DNA — `.fna` | `genome_sequence` | ✅ contigs w/ length/GC (+ `sequence` field available) |
| **Genotypic AMR / specialty genes** (new) | `sp_gene` | ✅ 574 hits for one genome; `property=Antibiotic Resistance` filterable |

**Implications:**

- **Protein FASTA for CD-HIT clustering can come from the API** — pull
  `genome_feature` (coords + `aa_sequence_md5`), then batch-fetch sequences from
  `feature_sequence` by md5. Removes Docker from the protein path. For thousands
  of genomes, dedupe md5s before fetching (identical proteins share an md5).
- **`sp_gene` is a new, relevant data source** — genotypic AMR determinants
  (resistance genes) to complement the phenotypic `genome_amr`. Not currently
  used by amRdata; worth considering for AMR modeling features.
- **Whole-assembly bulk files** (`.fna`) are the one case where **FTP is still
  the better route** — reconstructing multi-MB assemblies from per-contig API
  sequence fields is far heavier than downloading the flat file. FTP
  (`ftp.bvbrc.org/genomes/<id>/`) was **unreachable from the eval sandbox
  (network-blocked, HTTP 000), so verify from a real machine** — but it remains
  the documented, efficient path for sequence files. Keep FTP for `.fna`; use the
  API for metadata, features, protein sequences, and genotypic AMR.

**Net:** the API can replace the Docker/`p3-*` path for everything except bulk
assembly (`.fna`) downloads, while adding faceted subsampling and a genotypic-AMR
source the package doesn't currently have.

---

## 9. Per-species roster & the "adjustable row limit" request

Motivation: the `bvbrc` Python package caps `limit="max"` at 25,000, which is the
**server's** ceiling (see §4 — `limit(30000)` returns 25,000, and `start >= 25000`
is rejected). Goal: gather evidence to ask the package maintainer to stop
silently truncating at 25k.

**Critical framing (get this right in the request):** you cannot "just raise the
number." 25,000 is enforced by BV-BRC's API, not the package. The correct ask is
**automatic pagination** (keyset — §4) so a query that matches N > 25,000 rows
returns all N instead of a silent first 25k. Expose it as e.g. `limit="all"` or a
`max_rows` parameter that paginates under the hood. A bigger single `limit` value
will be clamped to 25k server-side and change nothing.

Roster produced by [`dev/bvbrc_species_roster.R`](../dev/bvbrc_species_roster.R)
→ [`dev/bvbrc_species_roster.csv`](../dev/bvbrc_species_roster.csv). Counts are
header-only (`Content-Range`); definitions:
`clean` = `genome_quality = Good`; `amr_rows` = `genome_amr` rows via
`eq(genome_name, <species>)`.

**Headline:** of 25 WHO(2024)/CDC(2019) bacterial priority species,
**14 exceed the 25k AMR-row cap** and **7 exceed 25k in genome count** (so even
the genome-metadata pull truncates). Worst case **E. coli: 7,388,629 AMR rows —
a capped single pull returns 0.3% of the data.** Other 7-figure species:
M. tuberculosis 2.26M, K. pneumoniae 2.01M, S. enterica 1.97M, S. pneumoniae
1.06M. Full table in the CSV.

**Caveat on `clean_with_amr`:** it counts genomes whose `genome` record has the
`antimicrobial_resistance` summary field populated — which is **sparsely filled**
and undercounts genomes that actually have `genome_amr` phenotype rows (e.g.
*C. jejuni*: 478,916 AMR rows but only 96 flagged genomes). For an accurate
"clean genomes with phenotype data," intersect Good `genome_id`s with the distinct
`genome_id`s present in `genome_amr` (heavier; via native Solr faceting). Treat the
`clean_with_amr` column as a floor, not a true count.

### Recommended minimal `genome` column set ("optimize / ditch the rest")

All verified present on the `genome` collection; use in `select(...)`:

`genome_id`, `assembly_accession` (NCBI/GCA), `genbank_accessions`,
`genome_quality`, `genome_status`, `checkm_completeness`, `checkm_contamination`,
`cds`, `genome_length`, `gc_content`, `host_name`, `isolation_country`,
`geographic_group`, `species`, `taxon_id`.

## References
- BV-BRC Data API: https://www.bv-brc.org/api/doc/
- `bvbrc` Python package: https://pypi.org/project/bvbrc/ ·
  https://github.com/abates20/bvbrc ·
  https://bvbrc.readthedocs.io/en/latest/
