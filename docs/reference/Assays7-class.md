# The `Assays7` class

An S7 implementation of the biocspec `Assays` interface. Holds an
ordered list of matrix- or array-like measurement objects that all agree
on `dim()[1:2]` (clauses `SE-INV-1`/`SE-INV-2`).

## Usage

``` r
Assays7(data = list(), assay_names = character(0))
```

## Arguments

- data:

  A list of matrix- or array-like objects, all sharing the same
  `dim()[1:2]`.

- assay_names:

  Optional character vector of names for the elements of `data`, the
  same length as `data`.

## Examples

``` r
m1 <- matrix(1:6, nrow = 2, ncol = 3)
m2 <- matrix(7:12, nrow = 2, ncol = 3)
a <- Assays7(data = list(counts = m1, logcounts = m2), assay_names = c("counts", "logcounts"))
dim(a)
#> [1] 2 3
length(a)
#> [1] 2
```
