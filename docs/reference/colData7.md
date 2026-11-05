# Get or set the per-sample (column) annotation table

Get or set the per-sample (column) annotation table

## Usage

``` r
colData7(x)

colData7(x) <- value
```

## Arguments

- x:

  A `SummarizedExperiment7`.

- value:

  Replacement `data.frame`, one row per sample; `rownames(value)` become
  `colnames(x)` when supplied.

## Examples

``` r
se <- SummarizedExperiment7(assays = list(counts = matrix(1:6, 2, 3)))
colData7(se) <- data.frame(group = c("a", "a", "b"), row.names = c("s1", "s2", "s3"))
colData7(se)
#>    group
#> s1     a
#> s2     a
#> s3     b
colnames(se)
#> [1] "s1" "s2" "s3"
```
