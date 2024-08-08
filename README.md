
# hdrm

<!-- badges: start -->
<!-- badges: end -->

R package for performing tests on high dimensional repeated measure data
for one group \[1\] or multiple groups \[2\].

## Installation

The current version can be installed with:

``` r
## install devtools package
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
# install package
devtools::install_github("Schnieboli/hdrm", dependencies = TRUE)
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo Schnieboli/hdrm@HEAD
#> cli     (3.6.2  -> 3.6.3 ) [CRAN]
#> ps      (1.7.6  -> 1.7.7 ) [CRAN]
#> withr   (3.0.0  -> 3.0.1 ) [CRAN]
#> pkgload (1.3.4  -> 1.4.0 ) [CRAN]
#> Rcpp    (1.0.12 -> 1.0.13) [CRAN]
#> Installing 5 packages: cli, ps, withr, pkgload, Rcpp
#> Installiere Pakete nach 'C:/Users/nilsh/AppData/Local/Temp/Rtmpo3OxRi/temp_libpathaf84ea61502'
#> (da 'lib' nicht spezifiziert)
#> Paket 'cli' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'ps' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'withr' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'pkgload' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> Paket 'Rcpp' erfolgreich ausgepackt und MD5 Summen abgeglichen
#> 
#> Die heruntergeladenen Binärpakete sind in 
#>  C:\Users\nilsh\AppData\Local\Temp\RtmpE7dQQX\downloaded_packages
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#>          checking for file 'C:\Users\nilsh\AppData\Local\Temp\RtmpE7dQQX\remotes1de427e45e49\Schnieboli-hdrm-d199db8/DESCRIPTION' ...     checking for file 'C:\Users\nilsh\AppData\Local\Temp\RtmpE7dQQX\remotes1de427e45e49\Schnieboli-hdrm-d199db8/DESCRIPTION' ...   ✔  checking for file 'C:\Users\nilsh\AppData\Local\Temp\RtmpE7dQQX\remotes1de427e45e49\Schnieboli-hdrm-d199db8/DESCRIPTION' (457ms)
#>       ─  preparing 'hdrm': (522ms)
#>    checking DESCRIPTION meta-information ...     checking DESCRIPTION meta-information ...   ✔  checking DESCRIPTION meta-information
#> ─  cleaning src
#>       ─  checking for LF line-endings in source and make files and shell scripts (514ms)
#>       ─  checking for empty or unneeded directories
#>       ─  building 'hdrm_0.0.0.9000.tar.gz'
#>      
#> 
#> Installiere Paket nach 'C:/Users/nilsh/AppData/Local/Temp/Rtmpo3OxRi/temp_libpathaf84ea61502'
#> (da 'lib' nicht spezifiziert)
```

## One group case

A one group test can be performed by using the function . Generic
functions exist for the classes , and .

    #> 
    #>           One Group Repeated Measure
    #>         
    #> Analysis of 20 individuals (0 removed) and 30 dimensions: 
    #> W = 0.156712  f = 20.55013  p.value = 0.39819 
    #> Null-Hypothesis: custom 
    #> Convergence parameter τ = 0.0486615
    #>   group    value sex subject variable region dimension
    #> 1  SCC+ 2.804180   W       1        1      1         1
    #> 2  SCC+ 2.285004   W       1        1      2         2
    #> 3  SCC+ 1.473094   W       1        1      3         3
    #> 4  SCC+ 1.548939   W       1        1      4         4
    #> 5  SCC+ 2.160346   W       1        1      5         5
    #> 6  SCC+ 2.143525   W       1        1      6         6
    #> 
    #>           One Group Repeated Measure
    #>         
    #> Analysis of 160 individuals (0 removed) and 40 dimensions: 
    #> W = 113.5808  f = 1.164473  p.value = 0 
    #> Null-Hypothesis: equal time profile 
    #> Convergence parameter τ = 0.8587573

## Multiple group case

<div id="refs" class="references csl-bib-body">

<div id="ref-Pauly2015-jo" class="csl-entry">

<span class="csl-left-margin">\[1\]
</span><span class="csl-right-inline">Pauly M, Ellenberger D and Brunner
E 2015 Analysis of high-dimensional one group repeated measures designs
*Statistics (Ber.)* **49** 1243–61</span>

</div>

<div id="ref-Sattler2018-az" class="csl-entry">

<span class="csl-left-margin">\[2\]
</span><span class="csl-right-inline">Sattler P and Pauly M 2018
Inference for high-dimensional split-plot-designs: A unified approach
for small to large numbers of factor levels *Electron. J. Stat.* **12**
2743–805</span>

</div>

</div>
