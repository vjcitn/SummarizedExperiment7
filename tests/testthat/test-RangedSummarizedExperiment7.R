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

test_that("rowRanges7<- requires matching length and does not promote/demote the class", {
  rse <- RangedSummarizedExperiment7(assays = m1, rowRanges = list(1, 2))
  rowRanges7(rse) <- list(10, 20)
  expect_identical(rse@row_ranges, list(10, 20))
  expect_error(rowRanges7(rse) <- list(1))
})
