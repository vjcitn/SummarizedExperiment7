# Get the names of the assays

Get the names of the assays

## Usage

``` r
assayNames7(x)
```

## Arguments

- x:

  A `SummarizedExperiment7`.

## Examples

``` r
se <- SummarizedExperiment7(assays = list(counts = matrix(1:6, 2, 3)))
assayNames7(se)
#> [1] "counts"
```
