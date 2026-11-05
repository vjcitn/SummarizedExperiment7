m1 <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))

test_that("RangedSummarizedExperiment7 forbids row_names and non-empty row_data (RSE-INV-3)", {
  rse <- RangedSummarizedExperiment7(
    assays = m1,
    rowRanges = list(g1 = 1:3, g2 = 4:6),
    colData = data.frame(group = c("a", "b", "c"), row.names = c("s1", "s2", "s3"))
  )
  expect_identical(dim(rse), c(2L, 3L))
  expect_null(rse@row_names)
  expect_identical(ncol(rse@row_data), 0L)
})

test_that("length(rowRanges) must equal nrow(x) (RSE-INV-2)", {
  expect_error(RangedSummarizedExperiment7(assays = m1, rowRanges = list(1)))
})

test_that("subsetting an RSE7 slices rowRanges in lockstep", {
  rse <- RangedSummarizedExperiment7(assays = m1, rowRanges = list(g1 = 1:3, g2 = 4:6))
  sub <- rse[1, 1:2]
  expect_identical(dim(sub), c(1L, 2L))
  expect_identical(sub@row_ranges, list(g1 = 1:3))
  expect_true(S7::S7_inherits(sub, RangedSummarizedExperiment7))
})

test_that("toRangedSummarizedExperiment7()/toSummarizedExperiment7() round-trip explicitly (RSE-COERCE-1/2)", {
  se <- SummarizedExperiment7(assays = m1, rowNames = c("g1", "g2"))
  rse <- toRangedSummarizedExperiment7(se, rowRanges = list(1, 2))
  expect_true(S7::S7_inherits(rse, RangedSummarizedExperiment7))
  back <- toSummarizedExperiment7(rse)
  expect_true(S7::S7_inherits(back, SummarizedExperiment7))
  expect_false(S7::S7_inherits(back, RangedSummarizedExperiment7))
})

test_that("toSummarizedExperiment7()/toRangedSummarizedExperiment7() reject non-S7 input instead of silently passing it through", {
  expect_error(toSummarizedExperiment7(1:3))
  expect_error(toSummarizedExperiment7(list(a = 1)))
  expect_error(toRangedSummarizedExperiment7(1:3, rowRanges = list(1, 2)))
})

test_that("cbind() unions rowRanges mcols() across objects and rejects mismatched range coordinates (SE-CBIND-1, conformance SE-CBIND-001/002)", {
  skip_if_not_installed("GenomicRanges")
  skip_if_not_installed("S4Vectors")
  skip_if_not_installed("IRanges")

  gr1 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(100, 500), width = c(50, 100)),
    strand = c("+", "-")
  )
  names(gr1) <- c("g1", "g2")
  S4Vectors::mcols(gr1) <- S4Vectors::DataFrame(symbol = c("A", "B"))
  rse1 <- RangedSummarizedExperiment7(
    assays = matrix(c(1, 3, 2, 4), 2, 2), rowRanges = gr1,
    colData = data.frame(trt = c("ctrl", "trt"), row.names = c("s1", "s2"))
  )

  gr2 <- gr1
  S4Vectors::mcols(gr2) <- S4Vectors::DataFrame(symbol = c("A", "B"), biotype = c("pc", "lnc"))
  rse2 <- RangedSummarizedExperiment7(
    assays = matrix(c(5, 6), 2, 1), rowRanges = gr2,
    colData = data.frame(trt = c("ctrl"), row.names = c("s3"))
  )

  bound <- cbind(rse1, rse2)
  expect_identical(dim(bound), c(2L, 3L))
  expect_identical(colnames(bound), c("s1", "s2", "s3"))
  expect_identical(as.character(rowData7(bound)$biotype), c("pc", "lnc"))

  gr3 <- GenomicRanges::GRanges(
    seqnames = c("chr1", "chr2"),
    ranges = IRanges::IRanges(start = c(999, 500), width = c(50, 100)),
    strand = c("+", "-")
  )
  names(gr3) <- c("g1", "g2")
  rse3 <- RangedSummarizedExperiment7(
    assays = matrix(c(5, 6), 2, 1), rowRanges = gr3,
    colData = data.frame(trt = c("ctrl"), row.names = c("s3"))
  )
  expect_error(cbind(rse1, rse3))
})

test_that("rowRanges7<- requires matching length and does not promote/demote the class", {
  rse <- RangedSummarizedExperiment7(assays = m1, rowRanges = list(1, 2))
  rowRanges7(rse) <- list(10, 20)
  expect_identical(rse@row_ranges, list(10, 20))
  expect_error(rowRanges7(rse) <- list(1))
})
