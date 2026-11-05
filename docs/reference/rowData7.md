# Get or set the per-feature (row) annotation table

Get or set the per-feature (row) annotation table

## Usage

``` r
rowData7(x)

rowData7(x) <- value
```

## Arguments

- x:

  A `SummarizedExperiment7`.

- value:

  Replacement `data.frame`, one row per feature.

## Examples

``` r
se <- SummarizedExperiment7(assays = list(counts = matrix(1:6, 2, 3)))
rowData7(se) <- data.frame(gene_biotype = c("protein_coding", "lncRNA"))
rowData7(se)
#>     gene_biotype
#> 1 protein_coding
#> 2         lncRNA
```
