# Implements the biocspec `Assays` interface (specification Section 4):
# a named or unnamed list of matrix/array-like objects that all agree on
# their first two dimensions.

#' The `Assays7` class
#'
#' An S7 implementation of the biocspec `Assays` interface. Holds an
#' ordered list of matrix- or array-like measurement objects that all
#' agree on `dim()[1:2]` (clauses `SE-INV-1`/`SE-INV-2`).
#'
#' @examples
#' m1 <- matrix(1:6, nrow = 2, ncol = 3)
#' m2 <- matrix(7:12, nrow = 2, ncol = 3)
#' a <- Assays7(data = list(counts = m1, logcounts = m2), assay_names = c("counts", "logcounts"))
#' dim(a)
#' length(a)
#'
#' @export
Assays7 <- S7::new_class(
  "Assays7",
  properties = list(
    data = S7::new_property(S7::class_list, default = quote(list())),
    assay_names = S7::new_property(S7::class_character | NULL, default = NULL)
  ),
  validator = function(self) {
    data <- self@data
    if (length(data) == 0) {
      return(NULL)
    }
    # Clause [SE-INV-1]: every element is matrix-like with ndim >= 2.
    for (a in data) {
      d <- dim(a)
      if (is.null(d) || length(d) < 2L) {
        return("all assays must be matrix-like objects with 2 (or more) dimensions")
      }
    }
    # Clause [SE-INV-2]: all assays share the same nrow and ncol.
    ref <- dim(data[[1]])[1:2]
    for (a in data) {
      if (!identical(dim(a)[1:2], ref)) {
        return("all assays must have the same nrow and ncol")
      }
    }
    if (!is.null(self@assay_names) && length(self@assay_names) != length(data)) {
      "@assay_names must have the same length as @data"
    }
  }
)

S7::method(dim, Assays7) <- function(x) {
  if (length(x@data) == 0) {
    c(0L, 0L)
  } else {
    dim(x@data[[1]])[1:2]
  }
}

S7::method(length, Assays7) <- function(x) length(x@data)

#' Slice an `Assays7` collection along both dimensions
#'
#' Applies the same row and column subscripts to every assay in lockstep,
#' mirroring the `slice2D()` algorithm of the biocspec `Assays` interface.
#'
#' @noRd
assays7_slice <- function(x, i, j) {
  sliced <- lapply(x@data, function(a) {
    extra <- rep(list(quote(expr = )), length(dim(a)) - 2L)
    do.call(`[`, c(list(a, i, j), extra, list(drop = FALSE)))
  })
  Assays7(data = sliced, assay_names = x@assay_names)
}

#' Row- or column-bind a list of `Assays7` objects
#'
#' Strict binding (clauses `SE-RBIND-4`/`SE-CBIND-3`): every object must
#' have the same number of assays, and if any object has names, all
#' objects must have the identical set of names.
#'
#' @noRd
assays7_bind <- function(objects, along_cols) {
  len0 <- length(objects[[1]])
  for (obj in objects) {
    if (length(obj) != len0) {
      stop2("the objects to bind must have the same number of assays")
    }
  }
  if (len0 == 0) {
    return(Assays7())
  }
  all_names <- lapply(objects, function(obj) obj@assay_names)
  has_names <- !all(vapply(all_names, is.null, logical(1)))
  bind_fun <- if (along_cols) cbind else rbind
  if (!has_names) {
    bound <- lapply(seq_len(len0), function(idx) {
      do.call(bind_fun, lapply(objects, function(obj) obj@data[[idx]]))
    })
    return(Assays7(data = bound, assay_names = NULL))
  }
  ref_names <- sort(objects[[1]]@assay_names)
  for (obj in objects) {
    if (is.null(obj@assay_names) || !identical(sort(obj@assay_names), ref_names)) {
      stop2("assays must have the same names()")
    }
  }
  unique_names <- objects[[1]]@assay_names
  bound <- lapply(unique_names, function(nm) {
    do.call(bind_fun, lapply(objects, function(obj) obj@data[[match(nm, obj@assay_names)]]))
  })
  Assays7(data = bound, assay_names = unique_names)
}

stop2 <- function(...) stop(..., call. = FALSE)
