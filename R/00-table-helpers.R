# Internal helpers treating a plain data.frame as the biocspec `Table`
# interface (Section 3.2). Kept out of the public API: callers interact
# with rowData()/colData(), never with these directly.

table_empty <- function(nrow = 0L) {
  as.data.frame(matrix(nrow = nrow, ncol = 0L))
}

check_table <- function(x, arg_name) {
  if (!is.data.frame(x)) {
    return(sprintf("@%s must be a data.frame (the 'Table' interface)", arg_name))
  }
  NULL
}

# Zero-column data.frames don't preserve row counts under base R's `[`/
# `rbind()` (a long-standing base R quirk), so row-slicing and row-binding
# go through these two helpers rather than calling `[`/`rbind()` directly.
table_slice_rows <- function(df, i) {
  if (ncol(df) == 0) {
    out <- table_empty(length(i))
    rownames(out) <- rownames(df)[i]
    out
  } else {
    df[i, , drop = FALSE]
  }
}

table_rbind <- function(dfs) {
  if (ncol(dfs[[1]]) == 0) {
    out <- table_empty(sum(vapply(dfs, nrow, integer(1))))
    rn <- unlist(lapply(dfs, rownames))
    if (!is.null(rn)) rownames(out) <- rn
    out
  } else {
    do.call(rbind, dfs)
  }
}
