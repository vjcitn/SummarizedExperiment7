m1 <- matrix(1:6, nrow = 2, ncol = 3, dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
m2 <- matrix(as.numeric(1:6) / 10, nrow = 2, ncol = 3, dimnames = dimnames(m1))

test_that("construction derives dimnames from the first assay (SE-CTOR-4)", {
  se <- SummarizedExperiment7(assays = list(counts = m1, norm = m2))
  expect_identical(dim(se), c(2L, 3L))
  expect_identical(rownames(se), c("g1", "g2"))
  expect_identical(colnames(se), c("s1", "s2", "s3"))
  expect_identical(assayNames7(se), c("counts", "norm"))
})

test_that("assay shape must agree across assays (SE-INV-2)", {
  bad <- matrix(1:4, nrow = 2)
  expect_error(SummarizedExperiment7(assays = list(m1, bad)))
})

test_that("assay shape must agree with rowData/colData (SE-INV-3/4)", {
  cd <- data.frame(group = c("a", "b", "c"), row.names = c("s1", "s2", "s3"))
  rd <- data.frame(chr = c("1", "2"), row.names = NULL)
  se <- SummarizedExperiment7(assays = m1, rowData = rd, colData = cd)
  expect_identical(nrow(se), 2L)
  expect_identical(ncol(se), 3L)
  expect_error(SummarizedExperiment7(assays = m1, colData = cd[1:2, , drop = FALSE]))
})

test_that("assay7()/assays7() accessors work", {
  se <- SummarizedExperiment7(assays = list(counts = m1, norm = m2))
  expect_identical(assay7(se, "counts"), m1)
  expect_identical(assay7(se), m1)
  expect_identical(assays7(se)$norm, m2)
})

test_that("2D subsetting slices rows/cols in lockstep (SE-SUB-1)", {
  se <- SummarizedExperiment7(assays = list(counts = m1), colData = data.frame(
    group = c("a", "b", "c"), row.names = c("s1", "s2", "s3")
  ))
  sub <- se[1, c("s1", "s3")]
  expect_identical(dim(sub), c(1L, 2L))
  expect_identical(rownames(sub), "g1")
  expect_identical(colnames(sub), c("s1", "s3"))
  expect_identical(assay7(sub), m1[1, c(1, 3), drop = FALSE])
})

test_that("rbind requires identical colnames, cbind requires identical rownames", {
  se <- SummarizedExperiment7(assays = list(counts = m1))
  top <- se[1, ]
  bottom <- se[2, ]
  bound <- rbind(top, bottom)
  expect_identical(dim(bound), c(2L, 3L))
  expect_identical(rownames(bound), c("g1", "g2"))

  left <- se[, 1:2]
  right <- se[, 3]
  bound2 <- cbind(left, right)
  expect_identical(dim(bound2), c(2L, 3L))
  expect_identical(colnames(bound2), c("s1", "s2", "s3"))

  mismatched <- se[, 1]
  colnames(mismatched) <- "other"
  expect_error(rbind(se, mismatched))
})

test_that("zero-assay containers are legal (SE-INV-5)", {
  se <- SummarizedExperiment7(
    rowData = data.frame(chr = c("1", "2")),
    colData = data.frame(group = c("a", "b"), row.names = c("s1", "s2"))
  )
  expect_identical(dim(se), c(2L, 2L))
  expect_identical(assays7(se), list())
})

test_that("rowData7<-/colData7<-/metadata7<- mutate in place", {
  se <- SummarizedExperiment7(assays = m1)
  rowData7(se) <- data.frame(chr = c("1", "2"))
  expect_identical(rowData7(se)$chr, c("1", "2"))
  metadata7(se) <- list(study = "demo")
  expect_identical(metadata7(se)$study, "demo")
})

test_that("rbind() unions colData columns across objects instead of keeping only the first (SE-RBIND-1/3, conformance SE-RBIND-001)", {
  se1 <- SummarizedExperiment7(
    assays = list(counts = matrix(c(1, 3, 2, 4), 2, 2)),
    rowNames = c("g1", "g2"),
    colData = data.frame(batch = c("b1", "b1"), donor = c("d1", "d2"), row.names = c("s1", "s2"))
  )
  se2 <- SummarizedExperiment7(
    assays = list(counts = matrix(c(5, 6), 1, 2)),
    rowNames = c("g3"),
    colData = data.frame(batch = c("b1", "b1"), tissue = c("lung", "lung"), row.names = c("s1", "s2"))
  )
  r <- rbind(se1, se2)
  cd <- colData7(r)
  expect_identical(cd$donor, c("d1", "d2"))
  expect_identical(cd$tissue, c("lung", "lung"))
})

test_that("rbind() still rejects conflicting duplicate colData columns (SE-RBIND-3, conformance SE-RBIND-002)", {
  se1 <- SummarizedExperiment7(
    assays = matrix(c(1, 2), 1, 2),
    colData = data.frame(batch = c("b1", "b1"), row.names = c("s1", "s2"))
  )
  se2 <- SummarizedExperiment7(
    assays = matrix(c(3, 4), 1, 2),
    colData = data.frame(batch = c("CONFLICT_1", "CONFLICT_2"), row.names = c("s1", "s2"))
  )
  expect_error(rbind(se1, se2))
})
