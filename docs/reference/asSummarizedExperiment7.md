# Coerce a Bioconductor S4 `SummarizedExperiment` to its S7 equivalent

Builds a
[`SummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/SummarizedExperiment7.md)
(or, if `x` carries `rowRanges`, a
[`RangedSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/RangedSummarizedExperiment7.md))
from a real Bioconductor
`SummarizedExperiment`/`RangedSummarizedExperiment` S4 object, by
reading off `assays()`, `rowData()`/`rowRanges()`, `colData()`, and
`metadata()`. This is a one-way, read-only snapshot: the result does not
stay linked to `x`, and there is no S7-to-S4 inverse in this package.
Note that
[`toSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md)/[`toRangedSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md)
are a different operation (moving between this package's own two S7
classes); pass an S4 object to *this* function, not to those.

## Usage

``` r
asSummarizedExperiment7(x)
```

## Arguments

- x:

  A `SummarizedExperiment` or `RangedSummarizedExperiment` S4 object, as
  defined by the Bioconductor `SummarizedExperiment` package.

## Examples

``` r
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  se4 <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = matrix(1:6, nrow = 2, ncol = 3,
      dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))),
    colData = S4Vectors::DataFrame(
      group = c("a", "a", "b"), row.names = c("s1", "s2", "s3")
    )
  )
  se7 <- asSummarizedExperiment7(se4)
  se7
}
#> <SummarizedExperiment7> 2 x 3
#> assays: counts 
#> rownames: g1, g2 
#> colnames: s1, s2, s3 
```
