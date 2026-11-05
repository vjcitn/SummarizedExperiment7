# Implements the biocspec `SummarizedExperiment` interface (specification
# Section 5). Internal slot names below echo Section 5.2 of
# biocspec/SummarizedExperiment-specification.md:
#   [[assays]], [[NAMES]] (-> row_names), [[elementMetadata]] (-> row_data),
#   [[colData]] (-> col_data), [[metadata]].
# As in the real R/Bioconductor implementation, nrow(x) is *defined* by
# nrow(row_data) and ncol(x) is *defined* by nrow(col_data); assay shape
# is checked against those, not the other way around (clauses
# `SE-INV-3`/`SE-INV-4`).

as_assays7 <- function(assays) {
  if (is.null(assays)) {
    return(NULL)
  }
  if (S7::S7_inherits(assays, Assays7)) {
    return(assays)
  }
  if (is.matrix(assays) || is.array(assays)) {
    return(Assays7(data = list(assays), assay_names = NULL))
  }
  if (is.list(assays)) {
    return(Assays7(data = unname(assays), assay_names = names(assays)))
  }
  stop2("`assays` must be NULL, a matrix/array, a list of matrices, or an Assays7 object")
}

se7_validator <- function(self) {
  msg <- check_table(self@row_data, "row_data")
  if (!is.null(msg)) return(msg)
  msg <- check_table(self@col_data, "col_data")
  if (!is.null(msg)) return(msg)

  nr <- nrow(self@row_data)
  nc <- nrow(self@col_data)

  if (!is.null(self@row_names) && length(self@row_names) != nr) {
    return("@row_names must have the same length as @row_data has rows")
  }
  if (!is.null(self@assays)) {
    ad <- dim(self@assays)
    # Clauses [SE-INV-3]/[SE-INV-4]: assay shape must agree with the
    # row/column metadata tables, not the other way around.
    if (ad[1] != nr) {
      return("nrow of assays must equal nrow of rowData")
    }
    if (ad[2] != nc) {
      return("ncol of assays must equal nrow of colData")
    }
  }
  NULL
}

se7_constructor <- function(
  assays = NULL,
  rowData = NULL,
  colData = NULL,
  metadata = list(),
  rowNames = NULL
) {
  assays7 <- as_assays7(assays)

  # Clause [SE-CTOR-4]: fall back to dimnames taken from the first assay
  # when rowData/colData/rowNames are not supplied.
  dn <- if (!is.null(assays7) && length(assays7@data) > 0) dimnames(assays7@data[[1]]) else NULL

  if (is.null(rowData)) {
    nr <- if (!is.null(assays7)) dim(assays7)[1] else 0L
    rowData <- table_empty(nr)
  }
  if (is.null(rowNames) && !is.null(dn) && !is.null(dn[[1]])) {
    rowNames <- dn[[1]]
  }
  if (is.null(colData)) {
    nc <- if (!is.null(assays7)) dim(assays7)[2] else 0L
    colData <- table_empty(nc)
    if (!is.null(dn) && !is.null(dn[[2]])) {
      rownames(colData) <- dn[[2]]
    }
  }

  S7::new_object(
    S7::S7_object(),
    assays = assays7,
    row_names = rowNames,
    row_data = rowData,
    col_data = colData,
    metadata = metadata
  )
}

#' The `SummarizedExperiment7` class
#'
#' An S7 implementation of the biocspec `SummarizedExperiment` interface:
#' a container coordinating one or more matrix-like assays with row
#' ("feature") and column ("sample") metadata tables.
#'
#' @param assays A matrix, list/named list of matrices, or an `Assays7`
#'   object. `NULL` for a zero-assay container.
#' @param rowData A `data.frame` of feature-level annotations, one row
#'   per feature.
#' @param colData A `data.frame` of sample-level annotations, one row
#'   per sample; `rownames(colData)` become `colnames(x)`.
#' @param metadata A named list of experiment-wide metadata.
#' @param rowNames Optional character vector of feature identifiers
#'   (`[[NAMES]]` in the specification); independent of `rowData`.
#'
#' @export
SummarizedExperiment7 <- S7::new_class(
  "SummarizedExperiment7",
  properties = list(
    assays = S7::new_property(Assays7 | NULL, default = NULL),
    row_names = S7::new_property(S7::class_character | NULL, default = NULL),
    row_data = S7::new_property(S7::class_data.frame, default = quote(table_empty(0L))),
    col_data = S7::new_property(S7::class_data.frame, default = quote(table_empty(0L))),
    metadata = S7::new_property(S7::class_list, default = quote(list()))
  ),
  validator = se7_validator,
  constructor = se7_constructor
)

## ---- Dimensions & identifiers (clauses SE-DIM-1/2) -----------------------

S7::method(dim, SummarizedExperiment7) <- function(x) {
  c(nrow(x@row_data), nrow(x@col_data))
}

S7::method(dimnames, SummarizedExperiment7) <- function(x) {
  list(x@row_names, rownames(x@col_data))
}

S7::method(`dimnames<-`, SummarizedExperiment7) <- function(x, value) {
  x@row_names <- value[[1]]
  rownames(x@col_data) <- value[[2]]
  x
}

S7::method(length, SummarizedExperiment7) <- function(x) nrow(x)
S7::method(names, SummarizedExperiment7) <- function(x) x@row_names
S7::method(`names<-`, SummarizedExperiment7) <- function(x, value) {
  x@row_names <- value
  x
}

## ---- Assay accessors & mutators (clause SE-ASSAY-1) -----------------------

#' @export
assays7 <- S7::new_generic("assays7", "x")
S7::method(assays7, SummarizedExperiment7) <- function(x) {
  if (is.null(x@assays)) list() else {
    nms <- x@assays@assay_names
    out <- x@assays@data
    if (!is.null(nms)) names(out) <- nms
    out
  }
}

#' @export
`assays7<-` <- S7::new_generic("assays7<-", "x")
S7::method(`assays7<-`, SummarizedExperiment7) <- function(x, value) {
  x@assays <- as_assays7(value)
  x
}

#' @export
assay7 <- S7::new_generic("assay7", "x")
S7::method(assay7, SummarizedExperiment7) <- function(x, i = 1L) {
  if (is.null(x@assays) || length(x@assays@data) == 0) {
    stop2("no assays present")
  }
  idx <- if (is.character(i)) match(i, x@assays@assay_names) else i
  if (is.na(idx)) stop2(sprintf("'%s' not in assayNames7(x)", i))
  x@assays@data[[idx]]
}

#' @export
`assay7<-` <- S7::new_generic("assay7<-", "x")
S7::method(`assay7<-`, SummarizedExperiment7) <- function(x, i = 1L, value) {
  a <- x@assays
  if (is.null(a)) a <- Assays7()
  idx <- if (is.character(i)) {
    found <- match(i, a@assay_names)
    if (is.na(found)) length(a@data) + 1L else found
  } else {
    i
  }
  data <- a@data
  data[[idx]] <- value
  nms <- a@assay_names
  if (is.character(i) && (length(nms) < idx || is.na(nms[idx]) || nms[idx] != i)) {
    if (length(nms) < idx) nms <- c(nms, rep(NA_character_, idx - length(nms)))
    nms[idx] <- i
  }
  x@assays <- Assays7(data = data, assay_names = nms)
  x
}

#' @export
assayNames7 <- S7::new_generic("assayNames7", "x")
S7::method(assayNames7, SummarizedExperiment7) <- function(x) {
  if (is.null(x@assays)) NULL else x@assays@assay_names
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## ---- Row, column, and experiment metadata (clauses SE-ROW/SE-COL/SE-META) --

#' @export
rowData7 <- S7::new_generic("rowData7", "x")
S7::method(rowData7, SummarizedExperiment7) <- function(x) x@row_data

#' @export
`rowData7<-` <- S7::new_generic("rowData7<-", "x")
S7::method(`rowData7<-`, SummarizedExperiment7) <- function(x, value) {
  x@row_data <- value
  x
}

#' @export
colData7 <- S7::new_generic("colData7", "x")
S7::method(colData7, SummarizedExperiment7) <- function(x) x@col_data

#' @export
`colData7<-` <- S7::new_generic("colData7<-", "x")
S7::method(`colData7<-`, SummarizedExperiment7) <- function(x, value) {
  rn <- rownames(x@col_data)
  x@col_data <- value
  if (is.null(rownames(value))) rownames(x@col_data) <- rn
  x
}

#' @export
metadata7 <- S7::new_generic("metadata7", "x")
S7::method(metadata7, SummarizedExperiment7) <- function(x) x@metadata

#' @export
`metadata7<-` <- S7::new_generic("metadata7<-", "x")
S7::method(`metadata7<-`, SummarizedExperiment7) <- function(x, value) {
  x@metadata <- value
  x
}

## ---- 2D subsetting (clauses SE-SUB-1..3) ----------------------------------

S7::method(`[`, SummarizedExperiment7) <- function(x, i, j, ..., drop = TRUE) {
  nr <- nrow(x)
  nc <- ncol(x)
  if (missing(i)) i <- seq_len(nr)
  if (missing(j)) j <- seq_len(nc)
  if (is.character(i)) i <- match(i, x@row_names)
  if (is.character(j)) j <- match(j, rownames(x@col_data))

  new_assays <- if (!is.null(x@assays)) assays7_slice(x@assays, i, j) else NULL
  new_row_data <- table_slice_rows(x@row_data, i)
  new_row_names <- if (!is.null(x@row_names)) x@row_names[i] else NULL
  new_col_data <- table_slice_rows(x@col_data, j)

  SummarizedExperiment7(
    assays = new_assays,
    rowData = new_row_data,
    colData = new_col_data,
    metadata = x@metadata,
    rowNames = new_row_names
  )
}

## ---- Combining (clauses SE-RBIND-*/SE-CBIND-*) ----------------------------

se7_rbind <- function(...) {
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
  SummarizedExperiment7(
    assays = new_assays,
    rowNames = unlist(lapply(objects, function(o) o@row_names)),
    rowData = table_rbind(lapply(objects, function(o) o@row_data)),
    colData = objects[[1]]@col_data,
    metadata = objects[[1]]@metadata
  )
}

se7_cbind <- function(...) {
  objects <- list(...)
  ref_rows <- objects[[1]]@row_names
  for (obj in objects) {
    if (!identical(obj@row_names, ref_rows)) {
      stop2("cbind: all objects must have identical rownames()")
    }
  }
  new_assays <- if (any(!vapply(objects, function(o) is.null(o@assays), logical(1)))) {
    assays7_bind(lapply(objects, function(o) o@assays %||% Assays7()), along_cols = TRUE)
  } else {
    NULL
  }
  SummarizedExperiment7(
    assays = new_assays,
    rowNames = objects[[1]]@row_names,
    rowData = objects[[1]]@row_data,
    colData = table_rbind(lapply(objects, function(o) o@col_data)),
    metadata = objects[[1]]@metadata
  )
}

S7::method(rbind, SummarizedExperiment7) <- function(..., deparse.level = 1) se7_rbind(...)
S7::method(cbind, SummarizedExperiment7) <- function(..., deparse.level = 1) se7_cbind(...)

## ---- print ------------------------------------------------------------

S7::method(print, SummarizedExperiment7) <- function(x, ...) {
  d <- dim(x)
  cat(sprintf(
    "<%s> %d x %d\n", S7::S7_class(x)@name, d[1], d[2]
  ))
  cat("assays:", paste(assayNames7(x) %||% character(), collapse = ", "), "\n")
  cat("rownames:", if (is.null(x@row_names)) "NULL" else paste(utils::head(x@row_names, 3), collapse = ", "), "\n")
  cat("colnames:", if (is.null(rownames(x@col_data))) "NULL" else paste(utils::head(rownames(x@col_data), 3), collapse = ", "), "\n")
  invisible(x)
}
