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

# TRUE iff `df`'s rownames were explicitly set, as opposed to being base R's
# automatic "1", "2", ... labels for a row.names-less data.frame. A
# zero-column table_empty() always carries the latter, and concatenating
# several objects' automatic labels (each restarting at "1") would produce
# spurious duplicates with no real meaning -- so only genuinely custom
# rownames (e.g. colData's, which double as colnames(x)) are preserved.
table_has_custom_rownames <- function(df) .row_names_info(df) >= 0L

table_rbind <- function(dfs) {
  if (ncol(dfs[[1]]) == 0) {
    out <- table_empty(sum(vapply(dfs, nrow, integer(1))))
    if (table_has_custom_rownames(dfs[[1]])) {
      rownames(out) <- unlist(lapply(dfs, rownames))
    }
    out
  } else {
    do.call(base::rbind, dfs)
  }
}

# Column-wise union of tables that describe the *same* rows (used by rbind()
# to merge colData across objects: the samples are identical, but each
# object may carry a different subset of annotation columns for them).
# A column present in more than one table must agree everywhere it overlaps
# (clause SE-RBIND-3); otherwise this errors.
table_merge_cols <- function(dfs) {
  nr <- nrow(dfs[[1]])
  out <- table_empty(nr)
  if (table_has_custom_rownames(dfs[[1]])) rownames(out) <- rownames(dfs[[1]])
  for (df in dfs) {
    if (ncol(df) == 0) next
    for (cn in names(df)) {
      if (cn %in% names(out)) {
        if (!isTRUE(all.equal(out[[cn]], df[[cn]], check.attributes = FALSE))) {
          stop2(sprintf("column '%s' is duplicated and the data do not match", cn))
        }
      } else {
        out[[cn]] <- df[[cn]]
      }
    }
  }
  out
}
