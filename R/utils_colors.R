# Colour palettes shared across the amR suite's plots (kept in sync with
# amRviz's R/utils_colors.R so reports and the dashboard read as one system).

# AMR phenotype palette (R/S/I plus full-word + lowercase variants so it works
# regardless of how the column is encoded). Susceptible is intentionally
# neutral grey so Resistant amber stands out as the signal of interest.
PHENOTYPE_COLORS <- c(
  R = "#d4872a", Resistant = "#d4872a", resistant = "#d4872a",
  S = "#5b8db8", Susceptible = "#5b8db8", susceptible = "#5b8db8",
  I = "#e6ab80", Intermediate = "#e6ab80", intermediate = "#e6ab80"
)

# Categorical palette for metadata plots (hosts, isolation sources, etc.).
# Use `meta_palette(n)` at call sites rather than referencing META_COLORS
# directly, so plots handle more categories than length(META_COLORS) without
# falling back to NA / gray70.
META_COLORS <- c(
  "#4a6b8a", "#87ceeb", "#e6ab80", "#8b6b7a",
  "#4e9a9a", "#c4a35a", "#9b7fba", "#5b8db8",
  "#d4735e", "#d4872a", "#b5b5b5"
)

#' Categorical palette of `n` muted colors
#'
#' For `n <= length(META_COLORS)` returns the first `n` curated colors verbatim;
#' above that, interpolates via `grDevices::colorRampPalette()` so callers get a
#' full palette regardless of how many unique values their data has.
#'
#' @param n Number of colors to return.
#' @return A character vector of `n` hex colors (empty when `n <= 0`).
#' @keywords internal
#' @noRd
meta_palette <- function(n = length(META_COLORS)) {
  if (n <= 0) {
    return(character(0))
  }
  if (n <= length(META_COLORS)) {
    return(META_COLORS[seq_len(n)])
  }
  grDevices::colorRampPalette(META_COLORS)(n)
}
