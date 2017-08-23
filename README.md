# accuMUlate demonstration 

This directory contains files used to re-analuse a mutation accumulation 
experiment using `accuMUlate`. The entire analysis is encapsulated in a
`Makefile` so running everting should only require you to type `make` in the
working directory. NOTE however, that a complete analysis will generate ~50Gb
of data

## Prerequesites

* [`accuMUlate`](https://github.com/dwinter/accuMUlate) and accuMUlate-tools.
* samtools 
* picard-tols
* gatk
* bedtools
* gnu parallel
* python

To generate a report describing these results you will also need  `R` and the
followign packages

* ggplot2
* reshape2
* rmarkdown

To  generate a PDF you will need a `TeX` environment and `pandoc`.

# Run everthing

```
$ make nproc=10 all
```
