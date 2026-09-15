clean_drug <- readr::read_tsv("data_raw/clean_drug.tsv")
drug_class <- readr::read_tsv("data_raw/drug_class.tsv")
drug_abbr <- readr::read_tsv("data_raw/drug_abbr.tsv")
class_abbr <- readr::read_tsv("data_raw/class_abbr.tsv")
cleaned_bvbrc_countries <- readr::read_tsv("data_raw/cleaned_bvbrc_countries.tsv")

usethis::use_data(
  clean_drug,
  drug_class,
  drug_abbr,
  class_abbr,
  cleaned_bvbrc_countries,
  overwrite = TRUE
)
