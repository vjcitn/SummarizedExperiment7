# SummarizedExperiment7: an experimental, agentically coded, biocspec-conformant S7 implementation

## What this package is

`SummarizedExperiment7` is an **experimental** reimplementation of the
`SummarizedExperiment`/`RangedSummarizedExperiment` data structures
using the `S7` object system instead of Bioconductor’s S4 classes.

Two things make this package unusual, and worth stating plainly before
using it:

- **It was written agentically.** The class contracts, accessors,
  validators, and tests were produced by an AI coding agent working from
  a specification document rather than hand-written line by line by a
  human. The code has been reviewed and exercised, but it has not had
  the years of production use that the Bioconductor S4
  `SummarizedExperiment` has had.
- **It targets a specification, not the S4 source.** Rather than porting
  Bioconductor’s S4 implementation, this package was built to conform to
  [biocspec](https://github.com/biocspec)’s language-agnostic
  specification for `SummarizedExperiment`
  (`SummarizedExperiment-specification.md`). That specification
  describes the *contract* – what a constructor, accessor, or combining
  operation must guarantee – independent of R, S4, or any particular
  language’s object system. biocspec also ships a formal conformance
  test suite (`summarized_experiment_conformance.yaml`); this package
  includes an adapted runner for it
  (`inst/conformance/run_conformance_se7.R`) covering the cases that are
  in scope for its API (see below).

Treat this package as a demonstration that a biocspec-conformant
`SummarizedExperiment` can be built outside of S4, not as a drop-in
replacement for the Bioconductor class.

## Deliberate differences from the S4 original

A handful of S4 behaviors are intentionally not replicated, because they
are implementation details rather than requirements of the biocspec
contract:

- No automatic class promotion/demotion when assigning to `rowRanges<-`;
  use
  [`toRangedSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md)/[`toSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md)
  instead.
- No `withDimnames=`/`checkDimnames=` constructor and accessor
  arguments.
- No `combineRows()`/`combineCols()` with fill values, no range-algebra
  delegation (`shift()`, `subsetByOverlaps()`, …), no `seqinfo()`
  pruning.

These are documented in the relevant `R/` source files and reflected in
the conformance runner’s “out of scope” list.

## Converting a real Bioconductor object: the `airway` example

The most common way to get data into `SummarizedExperiment7` from
existing Bioconductor workflows is
[`asSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/asSummarizedExperiment7.md),
which takes a real S4
`SummarizedExperiment`/`RangedSummarizedExperiment` and builds its S7
equivalent by reading off `assays()`, `rowRanges()`, `colData()`, and
`metadata()`. This is a one-way, read-only snapshot – the result is not
linked back to the original S4 object.

``` r

library(SummarizedExperiment7)
```

``` r

data(airway, package = "airway")
airway
```

    ## Loading required namespace: SummarizedExperiment

    ## class: RangedSummarizedExperiment 
    ## dim: 63677 8 
    ## metadata(1): ''
    ## assays(1): counts
    ## rownames(63677): ENSG00000000003 ENSG00000000005 ... ENSG00000273492
    ##   ENSG00000273493
    ## rowData names(10): gene_id gene_name ... seq_coord_system symbol
    ## colnames(8): SRR1039508 SRR1039509 ... SRR1039520 SRR1039521
    ## colData names(9): SampleName cell ... Sample BioSample

``` r

se7 <- asSummarizedExperiment7(airway)
se7
```

    ## <RangedSummarizedExperiment7> 63677 x 8
    ## assays: counts 
    ## rowRanges class: CompressedGRangesList 
    ## colnames: SRR1039508, SRR1039509, SRR1039512

``` r

class(se7)
```

    ## [1] "SummarizedExperiment7::RangedSummarizedExperiment7"
    ## [2] "SummarizedExperiment7::SummarizedExperiment7"      
    ## [3] "S7_object"

``` r

dim(se7)
```

    ## [1] 63677     8

The S7 object exposes the same core information through this package’s
own accessors:

``` r

assayNames7(se7)
```

    ## [1] "counts"

``` r

head(assay7(se7, "counts")[, 1:3])
```

    ##                 SRR1039508 SRR1039509 SRR1039512
    ## ENSG00000000003        679        448        873
    ## ENSG00000000005          0          0          0
    ## ENSG00000000419        467        515        621
    ## ENSG00000000457        260        211        263
    ## ENSG00000000460         60         55         40
    ## ENSG00000000938          0          0          2

``` r

head(colData7(se7))
```

    ##            SampleName    cell   dex albut        Run avgLength Experiment
    ## SRR1039508 GSM1275862  N61311 untrt untrt SRR1039508       126  SRX384345
    ## SRR1039509 GSM1275863  N61311   trt untrt SRR1039509       126  SRX384346
    ## SRR1039512 GSM1275866 N052611 untrt untrt SRR1039512       126  SRX384349
    ## SRR1039513 GSM1275867 N052611   trt untrt SRR1039513        87  SRX384350
    ## SRR1039516 GSM1275870 N080611 untrt untrt SRR1039516       120  SRX384353
    ## SRR1039517 GSM1275871 N080611   trt untrt SRR1039517       126  SRX384354
    ##               Sample    BioSample
    ## SRR1039508 SRS508568 SAMN02422669
    ## SRR1039509 SRS508567 SAMN02422675
    ## SRR1039512 SRS508571 SAMN02422678
    ## SRR1039513 SRS508572 SAMN02422670
    ## SRR1039516 SRS508575 SAMN02422682
    ## SRR1039517 SRS508576 SAMN02422673

``` r

rowRanges7(se7)
```

    ## GRangesList object of length 63677:
    ## $ENSG00000000003
    ## GRanges object with 17 ranges and 2 metadata columns:
    ##        seqnames            ranges strand |   exon_id       exon_name
    ##           <Rle>         <IRanges>  <Rle> | <integer>     <character>
    ##    [1]        X 99883667-99884983      - |    667145 ENSE00001459322
    ##    [2]        X 99885756-99885863      - |    667146 ENSE00000868868
    ##    [3]        X 99887482-99887565      - |    667147 ENSE00000401072
    ##    [4]        X 99887538-99887565      - |    667148 ENSE00001849132
    ##    [5]        X 99888402-99888536      - |    667149 ENSE00003554016
    ##    ...      ...               ...    ... .       ...             ...
    ##   [13]        X 99890555-99890743      - |    667156 ENSE00003512331
    ##   [14]        X 99891188-99891686      - |    667158 ENSE00001886883
    ##   [15]        X 99891605-99891803      - |    667159 ENSE00001855382
    ##   [16]        X 99891790-99892101      - |    667160 ENSE00001863395
    ##   [17]        X 99894942-99894988      - |    667161 ENSE00001828996
    ##   -------
    ##   seqinfo: 722 sequences (1 circular) from an unspecified genome
    ## 
    ## ...
    ## <63676 more elements>

Because `airway` is a `RangedSummarizedExperiment`,
[`asSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/asSummarizedExperiment7.md)
produced a `RangedSummarizedExperiment7`, whose per-feature identity
comes from
[`rowRanges7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/rowRanges7.md)
(a real `GRanges`) rather than an independent `rowData`/`rowNames` pair:

``` r

is(se7, "RangedSummarizedExperiment7")
```

    ## [1] FALSE

``` r

names(se7)[1:3]
```

    ## [1] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419"

Ordinary `SummarizedExperiment7` operations – subsetting,
[`rbind()`](https://rdrr.io/r/base/cbind.html),
[`cbind()`](https://rdrr.io/r/base/cbind.html) – work on the result
exactly as they would on an object built directly with
[`SummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/SummarizedExperiment7.md)/[`RangedSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/RangedSummarizedExperiment7.md):

``` r

se7[1:5, 1:2]
```

    ## <RangedSummarizedExperiment7> 5 x 2
    ## assays: counts 
    ## rowRanges class: CompressedGRangesList 
    ## colnames: SRR1039508, SRR1039509

## What `asSummarizedExperiment7()` is not

[`asSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/asSummarizedExperiment7.md)
converts *from* the S4 world *into* this package. It is a different
operation from
[`toSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md)/[`toRangedSummarizedExperiment7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/toRangedSummarizedExperiment7.md),
which move between this package’s own two S7 classes and will error if
given an S4 object – see their documentation for details. There is no
S7-to-S4 inverse in this package.

## Comparison to Bioconductor’s own S7 migration

Bioconductor’s `SummarizedExperiment` package has its own in-progress S7
migration, developed on the `s7-migration` branch of
[lawremi/SummarizedExperiment](https://github.com/lawremi/SummarizedExperiment/tree/s7-migration).
It is worth being explicit about how that effort differs from this
package, since both use S7 to represent the same data structure but for
different purposes.

That branch migrates the *actual* Bioconductor package to be built on S7
internally, while remaining a fully backward-compatible S4 object for
the roughly one thousand Bioconductor packages that depend on
`SummarizedExperiment`’s S4 API. Concretely: its core object is
constructed via
`new_object(S7_object(), NAMES=..., elementMetadata=..., colData=..., assays=..., metadata=...)`
and then wrapped into a real S4 object by `new_SummarizedExperiment()`;
every S7 generic (`rowData`, `colData`, `assays`, `assay`, `assayNames`,
and their setters) is also registered as a formal S4 generic via a
`setS4Generic()` helper, so existing `setMethod()`-based extensions in
dependent packages keep working untouched. It retains the complete
production feature set: `withDimnames=`/`checkDimnames=` arguments,
`GenomicRanges`/ `DelayedArray`/`HDF5Array`/`Seqinfo` integration,
`findOverlaps()`/ `nearest()`/`coverage()` methods, and the various
`makeSummarizedExperimentFrom*()` constructors.

`SummarizedExperiment7`, by contrast, is a small, standalone S7 class
hierarchy with no S4 inheritance and no S4/S7 bridging. It was built to
satisfy biocspec’s language-agnostic contract rather than to be a
drop-in replacement, and it deliberately omits the production features
listed above (see “Deliberate differences from the S4 original”). Its
accessors are namespaced
([`assay7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/assay7.md),
[`colData7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/colData7.md),
…) specifically so they never collide with the real package’s API.

|  | `SummarizedExperiment7` (this package) | `lawremi/SummarizedExperiment@s7-migration` |
|----|----|----|
| Purpose | Demonstrate a biocspec-conformant SE built from scratch in S7 | Migrate the real Bioconductor package’s internals to S7 |
| Relationship to S4 | None – pure S7, no S4 classes involved | Dual-natured: an S7 object wrapped as a real S4 object (`new_SummarizedExperiment()`) |
| Backward compatibility | Not applicable – new, namespaced API ([`assay7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/assay7.md), [`colData7()`](https://vjcitn.github.io/SummarizedExperiment7/reference/colData7.md), …) | Required – every S7 generic is also registered as an S4 generic (`setS4Generic()`) so dependent packages’ `setMethod()`s keep working |
| `withDimnames=`/`checkDimnames=` | Not implemented | Implemented |
| Range algebra (`shift()`, `subsetByOverlaps()`, `findOverlaps()`, `nearest()`, `coverage()`) | Not implemented | Implemented |
| `DelayedArray`/`HDF5Array`/`Seqinfo` integration | Not implemented | Implemented |
| `combineRows()`/`combineCols()` with fill | Not implemented | Not shown to be ported yet in the reviewed branch |
| [`rbind()`](https://rdrr.io/r/base/cbind.html) merges `colData` columns across objects | Yes (fixed this session; matches the real package’s `.cbind.DataFrame(args, colData, "colData")`) | Yes |
| [`cbind()`](https://rdrr.io/r/base/cbind.html) merges `rowRanges` `mcols()` across objects | Yes (fixed this session) | Yes |
| [`cbind()`](https://rdrr.io/r/base/cbind.html) rejects mismatched range coordinates, not just names/length | Yes (fixed this session) | Yes (`.compare(lapply(args, rowRanges), TRUE)`) |
| Code size (core class files) | ~3 files, a few hundred lines | ~2100 lines across `SummarizedExperiment-class.R`, `Assays-class.R`, `RangedSummarizedExperiment-class.R` alone |
| S7 version required | Any released S7 (CRAN 0.2.2 tested) | An unreleased S7 fork/branch – `R/acompat.R`/`R/Assays-class.R` call an `S4_contains()` helper and pass an S4 `classRepresentation` to `S7::new_class(parent=)`, neither of which CRAN S7 0.2.2 supports |
| biocspec conformance suite applied | Yes, 16/16 in-scope cases pass (`inst/conformance/run_conformance_se7.R`) | Not executable against this branch – see below (reviewed by reading source instead) |

Three of the real defects fixed in this package during development
([`rbind()`](https://rdrr.io/r/base/cbind.html) only keeping the first
object’s `colData` columns,
[`cbind()`](https://rdrr.io/r/base/cbind.html) not merging `rowRanges`
`mcols()`, and [`cbind()`](https://rdrr.io/r/base/cbind.html)’s
range-compatibility check only comparing names/length instead of actual
coordinates) turned out to match exactly how the canonical Bioconductor
S7 migration implements the same three operations – useful independent
validation that the fixes were in the right direction.

### Identifying and installing the matching S7 fork

`lawremi/SummarizedExperiment@s7-migration`’s most recent commit at the
time of writing (`f0689b1`, 2026-06-20 15:27:08 -0700) is titled “sync
with S7 branch API changes,” confirming it tracks an unreleased S7
branch rather than a CRAN release. That branch turned out to be
[RConsortium/S7@s4-extends-s7-compat-fixes](https://github.com/RConsortium/S7/tree/s4-extends-s7-compat-fixes).
(Checked first: lawremi’s own fork of S7, `lawremi/S7@local_fixes`,
which is older, 2024-era, and does not have the features `acompat.R`
needs.) The tip of that branch at the time of writing,
`c582779b4aa8cd1b50c4ca45e1e9b1083d342e04` (“Reduce S4 compatibility
test duplication,” 2026-06-25, version `0.2.2.9000`), was installed over
the previously installed CRAN S7 0.2.2 and does supply the
`S4_contains()` helper and S4-`classRepresentation`-as-`parent` support
that `s7-migration`’s `R/acompat.R` needs – but with it installed,
`lawremi/SummarizedExperiment@s7-migration` still fails to load, now at
`R/Assays-class.R`:
`Error in S4_contains(class): S4 classes can not extend abstract S7 classes.`
(`Assays` is declared `new_class(..., abstract = TRUE)` and then run
through `setS4Class()`/`S4_contains()`.)

That restriction was added to the S7 branch by commit `02a5348` (“Fix S4
bridge review findings,” 2026-06-24 09:14:37 -0400) and has held
unchanged through `c582779` (the tip). The commit that originally
introduced `S4_contains()`, `e75b66b` (2026-06-20 14:58:46 -0700), and
the next one after it, `abf31be` (“fix regression from introduction of
external generics,” 2026-06-20 15:26:01 -0700 – by Michael Lawrence,
just 67 seconds before `s7-migration`’s own sync commit, and almost
certainly the actual commit it was synced against), both predate that
restriction – but neither installs into a loadable package on its own
(`R CMD INSTALL` succeeds, but the freshly installed package fails
`.onLoad` with `no applicable method for` `` `@` ``
`` applied to an object of class "S7_S3_class"` ``), because the S7
branch was mid-way through a larger internal refactor (“massive
simplification” of its S3/S4 class representation, commit `67e88ae`) at
that exact point in its history. In other words: there is no commit of
`RConsortium/S7@s4-extends-s7-compat-fixes` that is *both* installable
and compatible with `s7-migration`’s current `Assays-class.R` – the
window where `S4_contains()` existed without the abstract-class
restriction was also a window where the package itself didn’t load. This
is a real, present incompatibility between the two branches at their
respective current tips, not a mistake in fork selection; S7’s own
installed copy was left at the loadable, restriction-bearing `c582779`,
since that is what `SummarizedExperiment7` itself was verified against
(below).

Installing `c582779` also surfaced and fixed two latent bugs in
`SummarizedExperiment7` itself: that commit’s revised method-dispatch
mechanism represents `rbind`/`cbind`/`` `[` `` as non-function
`S7_external_generic` sentinel objects in a package’s namespace (instead
of literally re-binding base’s own closures there, as CRAN S7 0.2.2
did), and two internal helpers (`assays7_bind()` in `R/01-assays.R`,
`table_rbind()` in `R/00-table-helpers.R`) had passed the bare symbol to
[`do.call()`](https://rdrr.io/r/base/do.call.html), which fetches the
binding directly and received the inert sentinel rather than a usable
function. Both now call
[`base::rbind()`](https://rdrr.io/r/base/cbind.html)/[`base::cbind()`](https://rdrr.io/r/base/cbind.html)/[`` base::`[` ``](https://rdrr.io/r/base/Extract.html)
explicitly. This also happens to be what resolves the “.Internal calls
in cbind/rbind” `R CMD check` WARNING noted elsewhere in this vignette’s
development history – the newer S7 fork’s sentinel-based dispatch no
longer exposes literal
[`.Internal()`](https://rdrr.io/r/base/Internal.html)-calling closures
under those names.

Because `s7-migration` cannot currently be loaded under any installable
S7 commit found, the biocspec conformance comparison in the table above
is based on reading its source (`R/SummarizedExperiment-class.R`,
`R/Assays-class.R`, `R/RangedSummarizedExperiment-class.R`), not on
executing its own test suite or the conformance runner against it.

As a check that the `base::`-qualified fix above is not itself a
regression, `SummarizedExperiment7` was also re-verified with the S7
fork uninstalled and the released CRAN S7 0.2.2 reinstalled in its
place: the same 62/62 testthat tests and 16/16 in-scope conformance
cases pass, and `R CMD check --as-cran --no-manual` returns to exactly
the pre-fork state of 1 WARNING (the inherent
[`.Internal()`](https://rdrr.io/r/base/Internal.html) call in
`cbind`/`rbind`, specific to CRAN S7 0.2.2’s dispatch mechanism) + 1
NOTE. The `base::` qualification is therefore safe under both S7
versions; the S7 fork was then reinstalled to match the rest of this
section.
