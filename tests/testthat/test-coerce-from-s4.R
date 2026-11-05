skip_if_not_installed("SummarizedExperiment")
skip_if_not_installed("S4Vectors")
skip_if_not_installed("GenomicRanges")
skip_if_not_installed("IRanges")

m <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))

test_that("asSummarizedExperiment7() rejects non-SummarizedExperiment input", {
  expect_error(asSummarizedExperiment7(1:3))
})

test_that("asSummarizedExperiment7() converts a plain SummarizedExperiment S4 object", {
  se4 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = m),
    rowData = S4Vectors::DataFrame(biotype = c("pc", "lnc"), row.names = c("g1", "g2")),
    colData = S4Vectors::DataFrame(group = c("a", "a", "b"), row.names = c("s1", "s2", "s3")),
    metadata = list(study = "demo")
  )
  se7 <- asSummarizedExperiment7(se4)
  expect_true(S7::S7_inherits(se7, SummarizedExperiment7))
  expect_false(S7::S7_inherits(se7, RangedSummarizedExperiment7))
  expect_identical(dim(se7), c(2L, 3L))
  expect_identical(rownames(se7), c("g1", "g2"))
  expect_identical(colnames(se7), c("s1", "s2", "s3"))
  expect_identical(assay7(se7, "counts"), m)
  expect_identical(rowData7(se7)$biotype, c("pc", "lnc"))
  expect_identical(metadata7(se7), list(study = "demo"))
})

test_that("asSummarizedExperiment7() converts a RangedSummarizedExperiment S4 object", {
  gr <- GenomicRanges::GRanges(
    c("chr1", "chr1"), IRanges::IRanges(c(1, 100), width = 10),
    biotype = c("pc", "lnc")
  )
  names(gr) <- c("g1", "g2")
  rse4 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = m),
    rowRanges = gr,
    colData = S4Vectors::DataFrame(group = c("a", "a", "b"), row.names = c("s1", "s2", "s3"))
  )
  rse7 <- asSummarizedExperiment7(rse4)
  expect_true(S7::S7_inherits(rse7, RangedSummarizedExperiment7))
  expect_identical(dim(rse7), c(2L, 3L))
  expect_identical(colnames(rse7), c("s1", "s2", "s3"))
  expect_s4_class(rowRanges7(rse7), "GRanges")
  expect_identical(names(rowRanges7(rse7)), c("g1", "g2"))
})

test_that("asSummarizedExperiment7() handles a zero-assay SummarizedExperiment S4 object", {
  se4 <- SummarizedExperiment::SummarizedExperiment(
    rowData = S4Vectors::DataFrame(x = 1:2),
    colData = S4Vectors::DataFrame(y = 1:3)
  )
  se7 <- asSummarizedExperiment7(se4)
  expect_identical(dim(se7), c(2L, 3L))
  expect_identical(assayNames7(se7), NULL)
})
