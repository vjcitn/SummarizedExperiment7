# Convert between `SummarizedExperiment7` and `RangedSummarizedExperiment7`

Explicit coercion functions, matching biocspec's
`toRangedSummarizedExperiment()` and `toSummarizedExperiment()` methods
(clauses `RSE-COERCE-1`/`RSE-COERCE-2`). Unlike the R/S4 original,
assigning to `rowRanges7<-` never implicitly changes an object's class;
these functions are the only way to move between the two classes. `x`
must already be one of this package's S7 classes — to convert a
Bioconductor S4 `SummarizedExperiment` or `RangedSummarizedExperiment`
object, use
[`asSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/asSummarizedExperiment7.md)
first (these functions error if `x` is not already S7).

## Usage

``` r
toRangedSummarizedExperiment7(x, rowRanges)

toSummarizedExperiment7(x)
```

## Arguments

- x:

  A `SummarizedExperiment7` (for `toRangedSummarizedExperiment7()`) or
  `RangedSummarizedExperiment7` (for `toSummarizedExperiment7()`).

- rowRanges:

  The `rowRanges7` to attach; its length must equal `nrow(x)`.

## Examples

``` r
se <- SummarizedExperiment7(assays = matrix(1:6, 2, 3), rowNames = c("g1", "g2"))
rse <- toRangedSummarizedExperiment7(se, rowRanges = list(1:10, 11:20))
rse
#> <RangedSummarizedExperiment7> 2 x 3
#> assays:  
#> rowRanges class: list 
#> colnames: 1, 2, 3 
back <- toSummarizedExperiment7(rse)
back
#> <SummarizedExperiment7> 2 x 3
#> assays:  
#> rownames: NULL 
#> colnames: 1, 2, 3 
```
