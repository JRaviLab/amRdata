# data_raw

Source tables for the package data in `data/*.rda`. After editing a table here,
rebuild the data with `source("data_raw/make_data.R")`.

## Drug tables

| File | Maps |
|---|---|
| `clean_drug.tsv` | raw BV-BRC antibiotic name to cleaned drug name (blank = not a drug, excluded) |
| `drug_class.tsv` | cleaned drug to drug class |
| `drug_abbr.tsv` | cleaned drug to abbreviation |
| `class_abbr.tsv` | drug class to abbreviation |

They are built from every antibiotic name in BV-BRC AMR records, with classes
suggested from BV-BRC (pharmacologic class, ATC), the WHO AWaRe list and
AntibioticDB (`ADB_all_compounds.csv`, DOI 10.1093/jac/dky208). The final tables
are only written after you review the proposals.

### Adding or updating drugs

Run from the package root (needs internet). The functions are in `R/drugclass.R`.

```r
devtools::load_all()
prepareDrugTables()      # 1. writes four *_review.tsv files
# 2. edit the review files (below)
finalizeDrugTables()     # 3. checks them, writes the final TSVs
source("data_raw/make_data.R")
```

**1. `prepareDrugTables()`** starts from the existing final TSVs (missing ones are
treated as empty, so it also works from scratch) and writes
`clean_drug_review.tsv`, `drug_class_review.tsv`, `drug_abbr_review.tsv` and
`class_abbr_review.tsv` in this folder. If review files already exist it
continues from them and keeps your edits; delete them to start over.

**2. Review.** Edit the review files by hand.

- A cell that says `TODO` needs your decision. The `suggestion` column beside it
  is the best guess (BV-BRC synonym, typo fix, combination, or a class from
  WHO / AntibioticDB / BV-BRC); copy it in or type your own. `note` says where a
  value came from.
- Values that are already filled in (exact name matches, classes that fit the
  existing scheme, beta-lactamase-inhibitor combinations, generated
  abbreviations) are proposals. Check them too.
- `cleaned_drug`: lower case. `_` joins words inside one drug name
  (`polymyxin_b`, `clavulanic_acid`); `-` only separates the drugs of a
  combination (`amoxicillin-clavulanic_acid`). A drug name must not be a class
  name (`macrolides` and `sulfa` are classes, not drugs). Leave it **blank** for
  values that are not a single antibiotic (`instrument`, `cephalosporin`).
- After deciding the names in `clean_drug_review.tsv`, run `prepareDrugTables()`
  again. It adds class and abbreviation rows for the drugs you just decided.
- `drug_class`: use an existing class where one fits. A new class is fine and
  gets a row in `class_abbr_review.tsv`.
- Codes (`drug_abbr`, `class_abbr`): upper case letters, digits and `-` only. No
  `_`: the resistance summary joins class codes with `_` and counts them.
  Generated codes are guesses; standard ones (CRO, TZB, ...) can't be derived.
- Rows for names you excluded can stay as they are.

**3. `finalizeDrugTables()`** checks the reviewed tables and writes the four final
TSVs only if every check passes; otherwise it lists what to fix and writes
nothing. On success it removes the review files. Checks:

- nothing left as `TODO` for a drug that is used
- keys and codes are unique
- names use only `a-z 0-9 _ -`, codes only `A-Z 0-9 -`
- every used drug has a class and an abbreviation, and every used class has one
- no drug is named after a class (`penicillin` and `tetracycline` are allowed:
  they are real drugs)

Rows for drugs that no cleaned name uses are kept with a warning, or dropped if
still `TODO`.
