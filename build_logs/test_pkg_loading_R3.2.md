```
R --version
```

```
R version 3.2.3 (2015-12-10) -- "Wooden Christmas-Tree"
Copyright (C) 2015 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under the terms of the
GNU General Public License versions 2 or 3.
For more information about these matters see
http://www.gnu.org/licenses/.
```

```
R
library(MutationalPatterns)
```

```
Loading required package: GenomicRanges
Loading required package: BiocGenerics
Loading required package: parallel
  
Attaching package: ‘BiocGenerics’
  
The following objects are masked from ‘package:parallel’:
  
    clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    clusterExport, clusterMap, parApply, parCapply, parLapply,
    parLapplyLB, parRapply, parSapply, parSapplyLB
  
The following objects are masked from ‘package:stats’:
  
    IQR, mad, xtabs
  
The following objects are masked from ‘package:base’:
  
    anyDuplicated, append, as.data.frame, as.vector, cbind, colnames,
    do.call, duplicated, eval, evalq, Filter, Find, get, grep, grepl,
    intersect, is.unsorted, lapply, lengths, Map, mapply, match, mget,
    order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    union, unique, unlist, unsplit

Loading required package: S4Vectors
Loading required package: stats4
Loading required package: IRanges
Loading required package: GenomeInfoDb
Loading required package: NMF
Loading required package: pkgmaker
Loading required package: registry

Attaching package: ‘pkgmaker’

The following object is masked from ‘package:S4Vectors’:

    new2

The following object is masked from ‘package:base’:

    isNamespaceLoaded

Loading required package: rngtools
Loading required package: cluster
NMF - BioConductor layer [OK] | Shared memory capabilities [NO: bigmemory] | Cores 15/16
  To enable shared memory capabilities, try: install.extras('
NMF
')

Attaching package: ‘NMF’

The following object is masked from ‘package:GenomicRanges’:

    compare

The following object is masked from ‘package:IRanges’:

    compare

The following objects are masked from ‘package:S4Vectors’:

    compare, nrun

Warning messages:
1: replacing previous import by ‘ggplot2::unit’ when loading ‘NMF’ 
2: replacing previous import by ‘ggplot2::arrow’ when loading ‘NMF’
```

```
devtools::session_info()
```

```
Session info -------------------------------------------------------------------
 setting  value
 version  R version 3.2.3 (2015-12-10)
 system   x86_64, linux-gnu
 ui       X11
 language (EN)
 collate  en_US.UTF-8
 tz       <NA>
 date     2016-09-02
  
Packages -----------------------------------------------------------------------
 package              * version  date       source
 AnnotationDbi          1.32.3   2016-05-13 Bioconductor
 Biobase              * 2.30.0   2016-05-11 Bioconductor
 BiocGenerics         * 0.16.1   2016-05-11 Bioconductor
 BiocParallel           1.4.3    2016-05-11 Bioconductor
 biomaRt                2.26.1   2016-05-13 Bioconductor
 Biostrings             2.38.4   2016-05-11 Bioconductor
 bitops                 1.0-6    2013-08-17 CRAN (R 3.2.3)
 BSgenome               1.38.0   2016-05-11 Bioconductor
 cluster              * 2.0.4    2016-04-18 CRAN (R 3.2.3)
 codetools              0.2-14   2015-07-15 CRAN (R 3.2.3)
 colorspace             1.2-6    2015-03-11 CRAN (R 3.2.3)
 DBI                    0.5      2016-08-11 CRAN (R 3.2.3)
 devtools               1.12.0   2016-06-24 CRAN (R 3.2.3)
 digest                 0.6.10   2016-08-02 CRAN (R 3.2.3)
 doParallel             1.0.10   2015-10-14 CRAN (R 3.2.3)
 foreach                1.4.3    2015-10-13 CRAN (R 3.2.3)
 futile.logger          1.4.3    2016-07-10 CRAN (R 3.2.3)
 futile.options         1.0.0    2010-04-06 CRAN (R 3.2.3)
 GenomeInfoDb         * 1.6.3    2016-05-11 Bioconductor  
 GenomicAlignments      1.6.3    2016-05-11 Bioconductor  
 GenomicFeatures        1.22.13  2016-05-13 Bioconductor  
 GenomicRanges        * 1.22.4   2016-08-30 Bioconductor  
 ggplot2                2.1.0    2016-03-01 CRAN (R 3.2.3)
 gridBase               0.4-7    2014-02-24 CRAN (R 3.2.3)
 gridExtra              2.2.1    2016-02-29 CRAN (R 3.2.3)
 gtable                 0.2.0    2016-02-26 CRAN (R 3.2.3)
 IRanges              * 2.4.8    2016-05-11 Bioconductor  
 iterators              1.0.8    2015-10-13 CRAN (R 3.2.3)
 lambda.r               1.1.9    2016-07-10 CRAN (R 3.2.3)
 magrittr               1.5      2014-11-22 CRAN (R 3.2.3)
 memoise                1.0.0    2016-01-29 CRAN (R 3.2.3)
 munsell                0.4.3    2016-02-13 CRAN (R 3.2.3)
 MutationalPatterns   * 0.2      2016-09-02 Bioconductor  
 NMF                  * 0.20.6   2015-05-26 CRAN (R 3.2.3)
 pkgmaker             * 0.22     2014-05-14 CRAN (R 3.2.3)
 plyr                   1.8.4    2016-06-08 CRAN (R 3.2.3)
 pracma                 1.9.3    2016-05-29 cran (@1.9.3) 
 quadprog               1.5-5    2013-04-17 CRAN (R 3.2.3)
 RColorBrewer           1.1-2    2014-12-07 CRAN (R 3.2.3)
 Rcpp                   0.12.6   2016-07-19 cran (@0.12.6)
 RCurl                  1.95-4.8 2016-03-01 CRAN (R 3.2.3)
 registry             * 0.3      2015-07-08 CRAN (R 3.2.3)
 reshape2               1.4.1    2014-12-06 CRAN (R 3.2.3)
 rngtools             * 1.2.4    2014-03-06 CRAN (R 3.2.3)
 Rsamtools              1.22.0   2016-07-10 Bioconductor  
 RSQLite                1.0.0    2014-10-25 CRAN (R 3.2.3)
 rstudioapi             0.6      2016-06-27 CRAN (R 3.2.3)
 rtracklayer            1.30.4   2016-07-10 Bioconductor  
 S4Vectors            * 0.8.11   2016-07-03 Bioconductor  
 scales                 0.4.0    2016-02-26 CRAN (R 3.2.3)
 stringi                1.1.1    2016-05-27 CRAN (R 3.2.3)
 stringr                1.1.0    2016-08-19 CRAN (R 3.2.3)
 SummarizedExperiment   1.0.2    2016-05-11 Bioconductor  
 VariantAnnotation      1.16.4   2016-08-19 Bioconductor  
 withr                  1.0.2    2016-06-20 CRAN (R 3.2.3)
 XML                    3.98-1.4 2016-03-01 CRAN (R 3.2.3)
 xtable                 1.8-2    2016-02-05 CRAN (R 3.2.3)
 XVector                0.10.0   2016-05-11 Bioconductor  
 zlibbioc               1.16.0   2016-05-11 Bioconductor
```

END
