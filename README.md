# BayesRep

**BayesRep** is an R package for the analysis of replication studies using Bayes
factors.

Pawel, S., Held, L. (2022). The sceptical Bayes factor for the assessment of
  replication success. *Journal of the Royal Statistical Society Series B:
  Statistical Methodology*. 84(3): 879-911.
  DOI:[10.1111/rssb.12491](https://doi.org/10.1111/rssb.12491)

## Installation

```r
## from GitHub
## install.packages("remotes") # requires remotes package
remotes::install_github(repo = "SamCH93/BayesRep")
```

## Usage

``` r
library("BayesRep")

## Original effect estimate and standard error
to <- 0.2
so <- 0.05

## Replication effect estimate and standard error 
tr <- 0.15
sr <- 0.04

## Compute and format sceptical Bayes factor
scepticalBF <- BFs(to = to, so = so, tr = tr, sr = sr)
formatBF(scepticalBF)

#> [1] "1/20"

## Compute and format replication Bayes factor
repBF <- BFr(to = to, so = so, tr = tr, sr = sr)
formatBF(repBF)

#> [1] "1/521"

## Plot posterior distribution of effect size
repPosterior(to = to, so = so, tr = tr, sr = sr)
```
![Plot of effect size posterior distribution based on original and replication data.](posterior.png)

<!-- png(filename = "posterior.png", width = 1200, height = 600, pointsize = 20); repPosterior(to = to, so = so, tr = tr, sr = sr, lwd = 1.5); dev.off() -->
