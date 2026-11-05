# Adapted conformance runner against biocspec's own formal conformance suite
# (bioconductor/biocspec conformance/summarized_experiment_conformance.yaml,
# suite version 0.1.0, 28 cases), translated to this package's S7 API.
#
# This is NOT a line-for-line port of biocspec's run_conformance_r.R (which
# builds real S4 SummarizedExperiment objects): each of the 16 in-scope
# cases below is hand-translated to call SummarizedExperiment7()/
# RangedSummarizedExperiment7() and this package's own accessors, using the
# same `given`/`when`/`then` values from the YAML case.
#
# 12 of the 28 cases are out of scope and not attempted: SE7 deliberately
# does not implement withDimnames (CTOR-004, DIMNAMES-001/002/003),
# checkDimnames (CTOR-005/006), a unified rowData+rowRanges constructor
# (CTOR-007), combineRows/combineCols with fill (COMBROWS-001/COMBCOLS-001),
# automatic class promotion/demotion (RSE-COERCE-001), range algebra
# delegation (RSE-RANGE-001), or seqinfo pruning (RSE-SEQINFO-001) -- see
# R/03-RangedSummarizedExperiment.R's header comment for the rationale.
#
# Last run result: 13/16 in-scope cases pass. 3 confirmed real failures:
#   SE-RBIND-001 - rbind() keeps only the first object's colData instead of
#     unioning colData columns across objects (se7_rbind() in
#     R/02-SummarizedExperiment.R does not merge colData at all).
#   SE-CBIND-001 - cbind() does not merge/union rowRanges mcols() across
#     objects (se7_cbind.RangedSummarizedExperiment7 keeps only the first
#     object's row_ranges/mcols).
#   SE-CBIND-002 - cbind()'s compatibility check on RangedSummarizedExperiment7
#     only compares rowRanges names/length, not actual range identity
#     (start/width/strand), so objects with the same names but different
#     coordinates are silently combined instead of rejected.

options(error = NULL)
devtools::load_all("/Users/vincentcarey/SAFE_CLAUDE/SummarizedExperiment7", quiet = TRUE)

results <- list()
record <- function(id, title, ok, note = "") {
  results[[length(results) + 1]] <<- list(id = id, title = title, ok = ok, note = note)
}
check <- function(id, title, expr) {
  msg <- NULL
  ok <- tryCatch({ expr(); TRUE }, error = function(e) { msg <<- conditionMessage(e); FALSE })
  record(id, title, ok, if (ok) "" else msg)
}
stopifnot_eq <- function(actual, expected, what) {
  if (!isTRUE(all.equal(actual, expected))) stop(sprintf("%s: expected %s, got %s", what, paste(expected, collapse=","), paste(actual, collapse=",")))
}

## SE-CTOR-001: zero-assay construction
check("SE-CTOR-001", "Empty SE construction with zero assays", function() {
  se <- SummarizedExperiment7()
  stopifnot_eq(dim(se), c(0L, 0L), "dim")
  stopifnot_eq(length(se), 0L, "length")
})

## SE-CTOR-002: single unnamed matrix, synthesized rowData/colData
check("SE-CTOR-002", "Single unnamed matrix auto promotion", function() {
  m <- matrix(c(1,2,3,4,5,6,7,8,9,10,11,12), nrow = 4, ncol = 3, byrow = TRUE)
  se <- SummarizedExperiment7(assays = m)
  stopifnot_eq(dim(se), c(4L, 3L), "dim")
  stopifnot_eq(length(se), 4L, "length")
  stopifnot_eq(length(assays7(se)), 1L, "assay count")
  stopifnot_eq(nrow(rowData7(se)), 4L, "rowData nrow")
  stopifnot_eq(ncol(rowData7(se)), 0L, "rowData ncol")
  stopifnot_eq(nrow(colData7(se)), 3L, "colData nrow")
  stopifnot_eq(ncol(colData7(se)), 0L, "colData ncol")
})

## SE-CTOR-003: dimnames fallback from first assay
check("SE-CTOR-003", "Dimnames fallback from first assay", function() {
  m <- matrix(c(10,30,50,20,40,60), nrow = 3, ncol = 2,
    dimnames = list(c("g1","g2","g3"), c("s1","s2")))
  se <- SummarizedExperiment7(assays = list(counts = m))
  stopifnot_eq(dim(se), c(3L, 2L), "dim")
  stopifnot_eq(rownames(se), c("g1","g2","g3"), "rownames")
  stopifnot_eq(names(se), c("g1","g2","g3"), "names")
  stopifnot_eq(colnames(se), c("s1","s2"), "colnames")
})

## SE-INV-001: mismatched assay dims rejected
check("SE-INV-001", "Validity rejects mismatched assay dims", function() {
  m1 <- matrix(1:4, 2, 2)
  m2 <- matrix(1:6, 2, 3)
  ok <- tryCatch({ SummarizedExperiment7(assays = list(a1 = m1, a2 = m2)); FALSE },
    error = function(e) TRUE)
  if (!ok) stop("expected an error but none was raised")
})

## SE-INV-002: assay nrow vs rowData nrow mismatch
check("SE-INV-002", "Validity rejects assay nrow vs rowData nrow mismatch", function() {
  m <- matrix(1:6, nrow = 3, ncol = 2)
  rd <- data.frame(gene_id = c("g1", "g2"))
  ok <- tryCatch({ SummarizedExperiment7(assays = m, rowData = rd); FALSE },
    error = function(e) TRUE)
  if (!ok) stop("expected an error but none was raised")
})

## SE-INV-003: assay ncol vs colData nrow mismatch
check("SE-INV-003", "Validity rejects assay ncol vs colData nrow mismatch", function() {
  m <- matrix(1:4, nrow = 2, ncol = 2)
  cd <- data.frame(sample_id = c("s1", "s2", "s3"))
  ok <- tryCatch({ SummarizedExperiment7(assays = m, colData = cd); FALSE },
    error = function(e) TRUE)
  if (!ok) stop("expected an error but none was raised")
})

## SE-ROW-001: rowData retrieval unaffected by row_names NA/duplicates
check("SE-ROW-001", "rowData unaffected by duplicate/NA row names", function() {
  m <- matrix(1:6, nrow = 3, ncol = 2)
  se <- SummarizedExperiment7(
    assays = m,
    rowNames = c("geneA", NA, "geneA"),
    rowData = data.frame(biotype = c("pc", "lnc", "pc"))
  )
  rd <- rowData7(se)
  stopifnot_eq(nrow(rd), 3L, "rowData nrow")
  stopifnot_eq(is.null(rownames(rd)) || all(rownames(rd) == as.character(1:3)), TRUE, "rowData rownames not gene-derived")
  stopifnot_eq(rd$biotype, c("pc", "lnc", "pc"), "biotype")
})

## SE-MUT-001: copy-on-change referential isolation
check("SE-MUT-001", "Copy-on-change referential isolation", function() {
  se <- SummarizedExperiment7(assays = matrix(1:4, 2, 2), rowNames = c("g1","g2"))
  se2 <- se
  names(se2) <- c("mut_g1", "mut_g2")
  assays7(se2) <- list(mut_counts = matrix(5:8, 2, 2))
  stopifnot_eq(names(se), c("g1", "g2"), "original row names")
  stopifnot_eq(assayNames7(se), NULL, "original assay names")
  stopifnot_eq(names(se2), c("mut_g1", "mut_g2"), "mutated row names")
  stopifnot_eq(assayNames7(se2), "mut_counts", "mutated assay names")
})

## SE-SUB-001: synchronized 2D subsetting with reordering
check("SE-SUB-001", "Synchronized 2D subsetting with reordering", function() {
  m <- matrix(c(11,12,13,21,22,23,31,32,33), nrow = 3, ncol = 3, byrow = TRUE)
  se <- SummarizedExperiment7(
    assays = list(counts = m), rowNames = c("g1","g2","g3"),
    colData = data.frame(cond = c("A","B","C"), row.names = c("s1","s2","s3")),
    rowData = data.frame(chr = c("chr1","chr2","chr3"))
  )
  sub <- se[c(3,1), c(2,3)]
  stopifnot_eq(dim(sub), c(2L, 2L), "dim")
  stopifnot_eq(rownames(sub), c("g3","g1"), "rownames")
  stopifnot_eq(colnames(sub), c("s2","s3"), "colnames")
  stopifnot_eq(unname(assay7(sub)), matrix(c(32,33,12,13), nrow=2, byrow=TRUE), "assay data")
  stopifnot_eq(rowData7(sub)$chr, c("chr3","chr1"), "rowData chr")
  stopifnot_eq(colData7(sub)$cond, c("B","C"), "colData cond")
})

## SE-SUB-002: single row/col never drops
check("SE-SUB-002", "Single row/column subsetting never drops", function() {
  m <- matrix(c(10,40,20,50,30,60), nrow = 2, ncol = 3)
  se <- SummarizedExperiment7(assays = list(counts = m), rowNames = c("g1","g2"),
    colData = data.frame(row.names = c("s1","s2","s3")))
  sub <- se["g2", 1]
  stopifnot_eq(dim(sub), c(1L, 1L), "dim")
  stopifnot_eq(rownames(sub), "g2", "rownames")
  stopifnot_eq(colnames(sub), "s1", "colnames")
})

## SE-SUB-003: subsetting by non-existent name errors
check("SE-SUB-003", "Subsetting by non-existent row name errors", function() {
  m <- matrix(1:4, 2, 2)
  se <- SummarizedExperiment7(assays = m, rowNames = c("g1","g2"))
  ok <- tryCatch({ se[c("g1","NON_EXISTENT_GENE"), ]; FALSE }, error = function(e) TRUE)
  # SE7's `[` uses match() which returns NA rather than erroring -- check actual behavior
  if (!ok) stop("no error raised (match() returned NA silently instead of erroring)")
})

## SE-SUB-004: subsetting to zero rows/cols yields valid 0-extent object
check("SE-SUB-004", "Subsetting to zero rows yields valid 0-extent SE", function() {
  m <- matrix(1:4, 2, 2)
  se <- SummarizedExperiment7(assays = m, rowNames = c("g1","g2"),
    colData = data.frame(row.names = c("s1","s2")))
  sub <- se[integer(0), ]
  stopifnot_eq(dim(sub), c(0L, 2L), "dim")
  stopifnot_eq(colnames(sub), c("s1","s2"), "colnames")
})

## SE-RBIND-001: rbind merges features, unions colData columns
check("SE-RBIND-001", "rbind merges features and unions colData columns", function() {
  se1 <- SummarizedExperiment7(
    assays = list(counts = matrix(c(1,3,2,4), 2, 2), norm = matrix(c(10,30,20,40), 2, 2)),
    rowNames = c("g1","g2"),
    colData = data.frame(batch = c("b1","b1"), donor = c("d1","d2"), row.names = c("s1","s2")),
    rowData = data.frame(symbol = c("A","B"))
  )
  se2 <- SummarizedExperiment7(
    assays = list(norm = matrix(c(50,60), 1, 2), counts = matrix(c(5,6), 1, 2)),
    rowNames = c("g3"),
    colData = data.frame(batch = c("b1","b1"), tissue = c("lung","lung"), row.names = c("s1","s2")),
    rowData = data.frame(symbol = c("C"))
  )
  r <- rbind(se1, se2)
  stopifnot_eq(dim(r), c(3L, 2L), "dim")
  stopifnot_eq(rownames(r), c("g1","g2","g3"), "rownames")
  stopifnot_eq(colnames(r), c("s1","s2"), "colnames")
  cd <- colData7(r)
  stopifnot_eq(cd$donor, c("d1","d2"), "colData donor (requires column union across objects)")
  stopifnot_eq(cd$tissue, c("lung","lung"), "colData tissue (requires column union across objects)")
})

## SE-RBIND-002: rbind rejects conflicting duplicate colData columns
check("SE-RBIND-002", "rbind rejects conflicting colData columns", function() {
  se1 <- SummarizedExperiment7(assays = matrix(c(1,2),1,2),
    colData = data.frame(batch = c("b1","b1"), row.names = c("s1","s2")))
  se2 <- SummarizedExperiment7(assays = matrix(c(3,4),1,2),
    colData = data.frame(batch = c("CONFLICT_1","CONFLICT_2"), row.names = c("s1","s2")))
  ok <- tryCatch({ rbind(se1, se2); FALSE }, error = function(e) TRUE)
  if (!ok) stop("expected a conflict error but none was raised (rbind does not merge/validate colData columns)")
})

## SE-CBIND-001: cbind combines samples with identical rowRanges, merges mcols
check("SE-CBIND-001", "cbind combines samples and merges rowRanges mcols", function() {
  skip <- !requireNamespace("GenomicRanges", quietly = TRUE) || !requireNamespace("S4Vectors", quietly = TRUE)
  if (skip) stop("SKIPPED: GenomicRanges/S4Vectors not installed")
  gr <- GenomicRanges::GRanges(
    seqnames = c("chr1","chr2"), ranges = IRanges::IRanges(start = c(100,500), width = c(50,100)),
    strand = c("+","-")
  )
  names(gr) <- c("g1","g2")
  S4Vectors::mcols(gr) <- S4Vectors::DataFrame(symbol = c("A","B"))
  rse1 <- RangedSummarizedExperiment7(assays = matrix(c(1,3,2,4), 2, 2), rowRanges = gr,
    colData = data.frame(trt = c("ctrl","trt"), row.names = c("s1","s2")))

  gr2 <- gr
  S4Vectors::mcols(gr2) <- S4Vectors::DataFrame(symbol = c("A","B"), biotype = c("pc","lnc"))
  rse2 <- RangedSummarizedExperiment7(assays = matrix(c(5,6), 2, 1), rowRanges = gr2,
    colData = data.frame(trt = c("ctrl"), row.names = c("s3")))

  r <- cbind(rse1, rse2)
  stopifnot_eq(dim(r), c(2L, 3L), "dim")
  stopifnot_eq(colnames(r), c("s1","s2","s3"), "colnames")
  stopifnot_eq(unname(assay7(r)), matrix(c(1,3,2,4,5,6), nrow=2), "assay data")
  rd <- rowData7(r)
  stopifnot_eq(as.character(rd$biotype), c("pc","lnc"), "merged rowRanges mcols biotype (requires mcols union across cbind inputs)")
})

## SE-CBIND-002: cbind rejects incompatible rowRanges
check("SE-CBIND-002", "cbind rejects RSE objects with incompatible rowRanges", function() {
  skip <- !requireNamespace("GenomicRanges", quietly = TRUE)
  if (skip) stop("SKIPPED: GenomicRanges not installed")
  gr1 <- GenomicRanges::GRanges(seqnames = c("chr1","chr2"),
    ranges = IRanges::IRanges(start = c(100,500), width = c(50,100)), strand = c("+","-"))
  rse1 <- RangedSummarizedExperiment7(assays = matrix(1:4, 2, 2), rowRanges = gr1,
    colData = data.frame(row.names = c("s1","s2")))
  gr2 <- GenomicRanges::GRanges(seqnames = c("chr1","chr2"),
    ranges = IRanges::IRanges(start = c(999,500), width = c(50,100)), strand = c("+","-"))
  rse2 <- RangedSummarizedExperiment7(assays = matrix(5:6, 2, 1), rowRanges = gr2,
    colData = data.frame(row.names = c("s3")))
  ok <- tryCatch({ cbind(rse1, rse2); FALSE }, error = function(e) TRUE)
  if (!ok) stop("expected error for mismatched ranges but none raised (cbind only checks names/length, not range identity)")
})

cat("\n==================== CONFORMANCE RESULTS (in-scope cases) ====================\n")
n_pass <- 0; n_fail <- 0
for (r in results) {
  status <- if (r$ok) "PASS" else "FAIL"
  if (r$ok) n_pass <- n_pass + 1 else n_fail <- n_fail + 1
  cat(sprintf("[%s] %-14s %s\n", status, r$id, r$title))
  if (!r$ok) cat(sprintf("        -> %s\n", r$note))
}
cat(sprintf("\n%d / %d passed (in-scope cases attempted)\n", n_pass, length(results)))

out_of_scope <- c(
  "SE-CTOR-004 (withDimnames not implemented)",
  "SE-CTOR-005 (checkDimnames constructor option not implemented)",
  "SE-CTOR-006 (checkDimnames constructor option not implemented)",
  "SE-CTOR-007 (no unified constructor accepting both rowData/rowRanges to conflict-check)",
  "SE-DIMNAMES-001 (withDimnames not implemented)",
  "SE-DIMNAMES-002 (withDimnames not implemented)",
  "SE-DIMNAMES-003 (withDimnames not implemented)",
  "SE-COMBROWS-001 (combineRows with fill not implemented)",
  "SE-COMBCOLS-001 (combineCols with fill not implemented)",
  "RSE-COERCE-001 (automatic class promotion/demotion on rowRanges<- not implemented)",
  "RSE-RANGE-001 (shift/subsetByOverlaps delegation not implemented)",
  "RSE-SEQINFO-001 (seqinfo pruning not implemented)"
)
cat(sprintf("\n%d cases out of scope (deliberate SE7 simplifications), not attempted:\n", length(out_of_scope)))
for (o in out_of_scope) cat(" -", o, "\n")
cat(sprintf("\nTotal suite cases: %d\n", length(results) + length(out_of_scope)))
