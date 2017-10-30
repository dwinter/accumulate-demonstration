# accuMUlate demonstration 

This directory contains files used to re-analuse a mutation accumulation 
experiment using `accuMUlate`. The entire analysis is encapsulated in a
`Makefile` so running everting should only require you to type `make` in the
working directory. NOTE however, that a complete analysis will generate ~50Gb
of data. 

By default all analyses will run on a single processor, which will take
some time. The number of processors to use can be set from the command line by
passing an value to `nproc`. For example, to run everything with 10 processors
you could type.

```
$ make nproc=10 all
```

## Prerequesites

Running this pipeline requires quite a few programs.You need [`accuMUlate`](https://github.com/dwinter/accuMUlate) 
and `denominate` (which is built at the same time as `accuMUlate`) to be
installed and on `$PATH`. Instructions for compiling these programs if given
at the [`accuMUlate` wiki](https://github.com/dwinter/accuMUlate/wiki).

In addition you will need: 

* samtools 
* bwa-mem
* picard-tols
* gatk
* bedtools
* gnu parallel
* python

To generate a report describing these results you will also need  `R` and the
following package:

* ggplot2
* reshape2
* rmarkdown
* bookdown

To  generate a PDF you will need a `TeX` environment and `pandoc`. If you don't
have these installed you can generate an html version of the report.

```sh
make athal_analysis.html
```


