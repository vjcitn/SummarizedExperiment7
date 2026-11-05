# Implements the biocspec `RangedSummarizedExperiment` interface
# (specification Section 6). `row_ranges` is duck-typed rather than
# hard-wired to Bioconductor's GenomicRanges: anything with a `length()`
# method is accepted, and a `names()`/`[`/`c()` method is used opportunistically
# when present. This keeps the package installable without a Bioconductor
# dependency while still working unmodified with a real GRanges/GRangesList
# if S4Vectors + GenomicRanges are available (see Suggests).
#
# Deliberate simplification vs. the R/S4 original (explicitly licensed by
# the task): there is no automatic class promotion/demotion triggered by
# assigning to `rowRanges7<-`. Use `toRangedSummarizedExperiment7()` and
# `toSummarizedExperiment7()` instead, matching biocspec's own
# `toRangedSummarizedExperiment()`/`toSummarizedExperiment()` method names
# (clauses `RSE-COERCE-1`/`RSE-COERCE-2`).

rse7_validator <- function(self) {
  msg <- se7_validator(self)
  if (!is.null(msg)) return(msg)
  if (!is.null(self@row_names)) {
    return("@row_names must be NULL on a RangedSummarizedExperiment7 (names come from @row_ranges)")
  }
  if (ncol(self@row_data) != 0) {
    return("@row_data must have zero columns on a RangedSummarizedExperiment7 (per-feature metadata lives on @row_ranges)")
  }
  if (length(self@row_ranges) != nrow(self@row_data)) {
    return("length(@row_ranges) must equal nrow(x)")
  }
  NULL
}

rse7_constructor <- function(
  assays = NULL,
  rowRanges,
  colData = NULL,
  metadata = list()
) {
  n <- length(rowRanges)
  parent <- SummarizedExperiment7(
    assays = assays, rowData = table_empty(n), colData = colData, metadata = metadata
  )
  # Clause RSE-INV-3: identity comes from @row_ranges, not @row_names, so
  # discard any rowNames the parent constructor inferred via the SE-CTOR-4
  # assay-dimnames fallback.
  parent@row_names <- NULL
  S7::new_object(parent, row_ranges = rowRanges)
}

#' The `RangedSummarizedExperiment7` class
#'
#' An S7 implementation of the biocspec `RangedSummarizedExperiment`
#' interface: a `SummarizedExperiment7` whose per-feature identity and
#' metadata are anchored to `rowRanges` rather than an independent
#' `rowData`/`rowNames` pair (clause `RSE-INV-3`).
#'
#' @param rowRanges Any object with a `length()` method, one element per
#'   feature (row); typically a `GenomicRanges::GRanges` or
#'   `GenomicRanges::GRangesList`, but duck-typed so a plain list or
#'   vector also works for the purposes of this demonstration package.
#' @inheritParams SummarizedExperiment7
#'
#' @examples
#' m <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
#' rse <- RangedSummarizedExperiment7(
#'   assays = m,
#'   rowRanges = list(g1 = 1:10, g2 = 11:20)
#' )
#' rse
#' dim(rse)
#' rowRanges7(rse)
#'
#' @export
RangedSummarizedExperiment7 <- S7::new_class(
  "RangedSummarizedExperiment7",
  parent = SummarizedExperiment7,
  properties = list(
    row_ranges = S7::new_property(S7::class_any)
  ),
  validator = rse7_validator,
  constructor = rse7_constructor
)

## ---- names/dimnames come from @row_ranges, not @row_names -----------------

S7::method(dimnames, RangedSummarizedExperiment7) <- function(x) {
  list(duck_names(x@row_ranges), rownames(x@col_data))
}

S7::method(names, RangedSummarizedExperiment7) <- function(x) duck_names(x@row_ranges)

duck_names <- function(r) {
  if (is.null(r)) return(NULL)
  tryCatch(names(r), error = function(e) NULL)
}

## ---- rowRanges accessor/mutator (clause RSE-DELEGATE-1) -------------------

#' Get or set the per-feature ranges of a `RangedSummarizedExperiment7`
#'
#' @param x A `RangedSummarizedExperiment7`.
#' @param value Replacement ranges; `length(value)` must equal `nrow(x)`.
#'   Reassigning never changes the object's class (see
#'   [toRangedSummarizedExperiment7()]).
#'
#' @examples
#' rse <- RangedSummarizedExperiment7(
#'   assays = matrix(1:6, 2, 3),
#'   rowRanges = list(1:10, 11:20)
#' )
#' rowRanges7(rse)
#' rowRanges7(rse) <- list(100:110, 200:220)
#' rowRanges7(rse)
#'
#' @export
rowRanges7 <- S7::new_generic("rowRanges7", "x")
S7::method(rowRanges7, RangedSummarizedExperiment7) <- function(x) x@row_ranges

#' @rdname rowRanges7
#' @export
`rowRanges7<-` <- S7::new_generic("rowRanges7<-", "x")
S7::method(`rowRanges7<-`, RangedSummarizedExperiment7) <- function(x, value) {
  if (length(value) != nrow(x)) {
    stop2("length(value) must equal nrow(x); use toRangedSummarizedExperiment7() to replace an incompatible rowRanges")
  }
  x@row_ranges <- value
  x
}

## rowData7 on a RSE7 surfaces per-feature metadata carried on @row_ranges
## when the duck-typed object supports mcols()/mcols<-() (S4Vectors::Vector
## objects, e.g. a real GRanges); otherwise it is the zero-column table
## required by [RSE-INV-3].
S7::method(rowData7, RangedSummarizedExperiment7) <- function(x) {
  if (has_mcols(x@row_ranges)) get_mcols(x@row_ranges) else x@row_data
}
S7::method(`rowData7<-`, RangedSummarizedExperiment7) <- function(x, value) {
  if (!has_mcols(x@row_ranges)) {
    stop2("rowData7<- requires @row_ranges to support mcols()/mcols<-() (e.g. a GenomicRanges::GRanges)")
  }
  x@row_ranges <- set_mcols(x@row_ranges, value)
  x
}

has_mcols <- function(r) {
  requireNamespace("S4Vectors", quietly = TRUE) && methods::is(r, "Vector")
}
get_mcols <- function(r) as.data.frame(S4Vectors::mcols(r))
set_mcols <- function(r, value) {
  S4Vectors::mcols(r) <- value
  r
}

## ---- 2D subsetting, overriding the row_names-based SummarizedExperiment7
## method so that @row_ranges is sliced in lockstep (clause RSE-INV-2) -----

S7::method(`[`, RangedSummarizedExperiment7) <- function(x, i, j, ..., drop = TRUE) {
  nr <- nrow(x)
  nc <- ncol(x)
  if (missing(i)) i <- seq_len(nr)
  if (missing(j)) j <- seq_len(nc)
  if (is.character(i)) i <- match(i, duck_names(x@row_ranges))
  if (is.character(j)) j <- match(j, rownames(x@col_data))

  new_assays <- if (!is.null(x@assays)) assays7_slice(x@assays, i, j) else NULL
  new_col_data <- table_slice_rows(x@col_data, j)

  RangedSummarizedExperiment7(
    assays = new_assays,
    rowRanges = x@row_ranges[i],
    colData = new_col_data,
    metadata = x@metadata
  )
}

## ---- Combining (clauses SE-RBIND-*/SE-CBIND-*, applied to row_ranges) ----

S7::method(rbind, RangedSummarizedExperiment7) <- function(..., deparse.level = 1) {
  objects <- list(...)
  ref_cols <- rownames(objects[[1]]@col_data)
  for (obj in objects) {
    if (!identical(rownames(obj@col_data), ref_cols)) {
      stop2("rbind: all objects must have identical colnames()")
    }
  }
  new_assays <- if (any(!vapply(objects, function(o) is.null(o@assays), logical(1)))) {
    assays7_bind(lapply(objects, function(o) o@assays %||% Assays7()), along_cols = FALSE)
  } else {
    NULL
  }
  combined_ranges <- do.call(c, lapply(objects, function(o) o@row_ranges))
  RangedSummarizedExperiment7(
    assays = new_assays,
    rowRanges = combined_ranges,
    colData = table_merge_cols(lapply(objects, function(o) o@col_data)),
    metadata = objects[[1]]@metadata
  )
}

# Compares the ranges themselves (seqnames/start/width/strand via
# GenomicRanges::granges(), which drops mcols), not just names/length, so
# that cbind() rejects objects whose coordinates differ even when their
# names and length happen to match (clause SE-CBIND-1). Falls back to
# all.equal() for duck-typed non-GRanges ranges (e.g. plain lists).
ranges_identical <- function(a, b) {
  if (!identical(duck_names(a), duck_names(b))) return(FALSE)
  if (length(a) != length(b)) return(FALSE)
  if (requireNamespace("GenomicRanges", quietly = TRUE) &&
      methods::is(a, "GenomicRanges") && methods::is(b, "GenomicRanges")) {
    return(identical(
      as.character(GenomicRanges::granges(a)),
      as.character(GenomicRanges::granges(b))
    ))
  }
  isTRUE(all.equal(a, b))
}

S7::method(cbind, RangedSummarizedExperiment7) <- function(..., deparse.level = 1) {
  objects <- list(...)
  ref_ranges <- objects[[1]]@row_ranges
  for (obj in objects) {
    if (!ranges_identical(obj@row_ranges, ref_ranges)) {
      stop2("cbind: all objects must have identical rowRanges()")
    }
  }
  new_assays <- if (any(!vapply(objects, function(o) is.null(o@assays), logical(1)))) {
    assays7_bind(lapply(objects, function(o) o@assays %||% Assays7()), along_cols = TRUE)
  } else {
    NULL
  }
  merged_ranges <- ref_ranges
  if (has_mcols(ref_ranges)) {
    merged_mcols <- table_merge_cols(lapply(objects, function(o) get_mcols(o@row_ranges)))
    merged_ranges <- set_mcols(ref_ranges, merged_mcols)
  }
  RangedSummarizedExperiment7(
    assays = new_assays,
    rowRanges = merged_ranges,
    colData = table_rbind(lapply(objects, function(o) o@col_data)),
    metadata = objects[[1]]@metadata
  )
}

## ---- Promotion / demotion (clauses RSE-COERCE-1/2) ------------------------

#' Convert between `SummarizedExperiment7` and `RangedSummarizedExperiment7`
#'
#' Explicit coercion functions, matching biocspec's `toRangedSummarizedExperiment()`
#' and `toSummarizedExperiment()` methods (clauses `RSE-COERCE-1`/`RSE-COERCE-2`).
#' Unlike the R/S4 original, assigning to `rowRanges7<-` never implicitly
#' changes an object's class; these functions are the only way to move
#' between the two classes. `x` must already be one of this package's S7
#' classes — to convert a Bioconductor S4 `SummarizedExperiment` or
#' `RangedSummarizedExperiment` object, use [asSummarizedExperiment7()]
#' first (these functions error if `x` is not already S7).
#'
#' @param x A `SummarizedExperiment7` (for `toRangedSummarizedExperiment7()`)
#'   or `RangedSummarizedExperiment7` (for `toSummarizedExperiment7()`).
#' @param rowRanges The `rowRanges7` to attach; its length must equal
#'   `nrow(x)`.
#'
#' @examples
#' se <- SummarizedExperiment7(assays = matrix(1:6, 2, 3), rowNames = c("g1", "g2"))
#' rse <- toRangedSummarizedExperiment7(se, rowRanges = list(1:10, 11:20))
#' rse
#' back <- toSummarizedExperiment7(rse)
#' back
#'
#' @export
toRangedSummarizedExperiment7 <- function(x, rowRanges) {
  if (!S7::S7_inherits(x, SummarizedExperiment7)) {
    stop2(
      "`x` must be a SummarizedExperiment7 or RangedSummarizedExperiment7; ",
      "for a Bioconductor S4 object, convert it first with asSummarizedExperiment7()"
    )
  }
  if (S7::S7_inherits(x, RangedSummarizedExperiment7)) return(x)
  if (length(rowRanges) != nrow(x)) {
    stop2("length(rowRanges) must equal nrow(x)")
  }
  RangedSummarizedExperiment7(
    assays = x@assays,
    rowRanges = rowRanges,
    colData = x@col_data,
    metadata = x@metadata
  )
}

#' @rdname toRangedSummarizedExperiment7
#' @export
toSummarizedExperiment7 <- function(x) {
  if (!S7::S7_inherits(x, SummarizedExperiment7)) {
    stop2(
      "`x` must be a SummarizedExperiment7 or RangedSummarizedExperiment7; ",
      "for a Bioconductor S4 object, convert it first with asSummarizedExperiment7()"
    )
  }
  if (!S7::S7_inherits(x, RangedSummarizedExperiment7)) return(x)
  SummarizedExperiment7(
    assays = x@assays,
    rowData = table_empty(nrow(x)),
    rowNames = duck_names(x@row_ranges),
    colData = x@col_data,
    metadata = x@metadata
  )
}

## ---- print ------------------------------------------------------------

S7::method(print, RangedSummarizedExperiment7) <- function(x, ...) {
  d <- dim(x)
  cat(sprintf("<%s> %d x %d\n", S7::S7_class(x)@name, d[1], d[2]))
  cat("assays:", paste(assayNames7(x) %||% character(), collapse = ", "), "\n")
  cat("rowRanges class:", paste(class(x@row_ranges), collapse = "/"), "\n")
  cat("colnames:", if (is.null(rownames(x@col_data))) "NULL" else paste(utils::head(rownames(x@col_data), 3), collapse = ", "), "\n")
  invisible(x)
}
