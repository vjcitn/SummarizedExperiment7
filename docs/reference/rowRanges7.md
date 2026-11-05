# Get or set the per-feature ranges of a `RangedSummarizedExperiment7`

Get or set the per-feature ranges of a `RangedSummarizedExperiment7`

## Usage

``` r
rowRanges7(x)

rowRanges7(x) <- value
```

## Arguments

- x:

  A `RangedSummarizedExperiment7`.

- value:

  Replacement ranges; `length(value)` must equal `nrow(x)`. Reassigning
  never changes the object's class (see
  [`toRangedSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md)).

## Examples

``` r
rse <- RangedSummarizedExperiment7(
  assays = matrix(1:6, 2, 3),
  rowRanges = list(1:10, 11:20)
)
rowRanges7(rse)
#> [[1]]
#>  [1]  1  2  3  4  5  6  7  8  9 10
#> 
#> [[2]]
#>  [1] 11 12 13 14 15 16 17 18 19 20
#> 
rowRanges7(rse) <- list(100:110, 200:220)
rowRanges7(rse)
#> [[1]]
#>  [1] 100 101 102 103 104 105 106 107 108 109 110
#> 
#> [[2]]
#>  [1] 200 201 202 203 204 205 206 207 208 209 210 211 212 213 214 215 216 217 218
#> [20] 219 220
#> 
```
