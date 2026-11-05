# SummarizedExperiment7

An [S7](https://rconsortium.github.io/S7/) implementation of the
`SummarizedExperiment` and `RangedSummarizedExperiment` data structures,
built from the language-agnostic
[biocspec](https://github.com/jkanche/biocspec) specification rather
than from Bioconductor’s S4 source. It is conformant with that
specification without necessarily replicating every implementation
detail of the S4 original — notably, there is no automatic class
promotion/demotion on `rowRanges7<-`; use
[`toRangedSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md)
and
[`toSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md)
instead.

## Installation

``` r

# install.packages("remotes")
remotes::install_github("vjcitn/SummarizedExperiment7")
```

## Usage

``` r

library(SummarizedExperiment7)

m <- matrix(1:6, nrow = 2, ncol = 3,
            dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))

se <- SummarizedExperiment7(
  assays = list(counts = m),
  colData = data.frame(group = c("a", "a", "b"), row.names = c("s1", "s2", "s3"))
)

se
dim(se)
assay7(se, "counts")
colData7(se)
```

## Testing

``` r

devtools::test()
```
