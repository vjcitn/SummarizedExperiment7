# Get or set a single assay by position or name

Get or set a single assay by position or name

## Usage

``` r
assay7(x, i = 1L)

assay7(x, i = 1L) <- value
```

## Arguments

- x:

  A `SummarizedExperiment7`.

- i:

  Assay position (integer) or name (character); defaults to `1L`, the
  first assay.

- value:

  Replacement matrix for the `i`-th (or named) assay.

## Examples

``` r
se <- SummarizedExperiment7(assays = list(counts = matrix(1:6, 2, 3)))
assay7(se)
#>      [,1] [,2] [,3]
#> [1,]    1    3    5
#> [2,]    2    4    6
assay7(se, "counts")
#>      [,1] [,2] [,3]
#> [1,]    1    3    5
#> [2,]    2    4    6
assay7(se, "logcounts") <- log1p(assay7(se, "counts"))
assayNames7(se)
#> [1] "counts"    "logcounts"
```
