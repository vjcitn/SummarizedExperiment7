# Bridges the Bioconductor S4 SummarizedExperiment/RangedSummarizedExperiment
# classes into their S7 equivalents here. Guarded on the Bioconductor
# 'SummarizedExperiment' package being installed (it is Suggests-only, not
# a hard dependency of this package).

#' Coerce a Bioconductor S4 `SummarizedExperiment` to its S7 equivalent
#'
#' Builds a [SummarizedExperiment7()] (or, if `x` carries `rowRanges`, a
#' [RangedSummarizedExperiment7()]) from a real Bioconductor
#' `SummarizedExperiment`/`RangedSummarizedExperiment` S4 object, by reading
#' off `assays()`, `rowData()`/`rowRanges()`, `colData()`, and `metadata()`.
#' This is a one-way, read-only snapshot: the result does not stay linked to
#' `x`, and there is no S7-to-S4 inverse in this package.
#'
#' @param x A `SummarizedExperiment` or `RangedSummarizedExperiment` S4
#'   object, as defined by the Bioconductor `SummarizedExperiment` package.
#'
#' @examples
#' if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
#'   se4 <- SummarizedExperiment::SummarizedExperiment(
#'     assays = list(counts = matrix(1:6, nrow = 2, ncol = 3,
#'       dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))),
#'     colData = S4Vectors::DataFrame(
#'       group = c("a", "a", "b"), row.names = c("s1", "s2", "s3")
#'     )
#'   )
#'   se7 <- asSummarizedExperiment7(se4)
#'   se7
#' }
#'
#' @export
asSummarizedExperiment7 <- function(x) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop2("asSummarizedExperiment7() requires the Bioconductor 'SummarizedExperiment' package")
  }
  if (!methods::is(x, "SummarizedExperiment")) {
    stop2("`x` must be a SummarizedExperiment (or RangedSummarizedExperiment) S4 object")
  }

  assay_list <- as.list(SummarizedExperiment::assays(x))
  if (length(assay_list) == 0) assay_list <- NULL

  col_data <- as.data.frame(SummarizedExperiment::colData(x))
  meta <- as.list(S4Vectors::metadata(x))

  if (methods::is(x, "RangedSummarizedExperiment")) {
    RangedSummarizedExperiment7(
      assays = assay_list,
      rowRanges = SummarizedExperiment::rowRanges(x),
      colData = col_data,
      metadata = meta
    )
  } else {
    SummarizedExperiment7(
      assays = assay_list,
      rowData = as.data.frame(SummarizedExperiment::rowData(x)),
      colData = col_data,
      metadata = meta,
      rowNames = rownames(x)
    )
  }
}
