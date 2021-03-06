# deplink
### an R package to compare the genetic/epigenetic features between cancer cell lines with different dependencies of a gene set (signature)


[![](https://img.shields.io/badge/release%20version-0.99.0-green.svg)](https://github.com/seanchen607/deplink)
[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

‘deplink’ compares the genetic/epigenetic features between cancer cell lines with different dependencies of a gene set (signature).

Data source: [DepMap](https://depmap.org/portal/) *(release 2019q4)* and [CCLE](https://portals.broadinstitute.org/ccle)

For details, please see [Tutorial](https://seanchen607.github.io/deplink.html).

<!--
<object data="docs/dep_9-1-1_score_CancerType.TCGA.dotplot.pdf" type="application/pdf" width="700px" height="200px">
    <embed src="docs/dep_9-1-1_score_CancerType.TCGA.dotplot.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="docs/dep_9-1-1_score_CancerType.TCGA.dotplot.pdf">Download PDF</a>.</p>
    </embed>
</object>

<embed src="https://drive.google.com/viewerng/viewer?embedded=true&url=https://github.com/seanchen607/deplink/raw/master/docs/dep_9-1-1_score_CancerType.TCGA.dotplot.pdf" type='application/pdf'>
-->

<!--
<a href="https://doi.org/10.1111/1755-0998.13023"><img src="docs/schematic.png" height="100" align="center" /></a>
-->

## :gear: Install deplink in R (>= 3.5.0)

	library(devtools)
	install_github("seanchen607/deplink")
	
## :hourglass_flowing_sand: Load deplink and data libraries

	library(deplink)
	source(system.file("script", "load_libs.R", package = "deplink"))

## :dna: Usage 

For example, deplink compares the genetic/epigenetic features between cancer cell lines with highest and lowest dependencies of "9-1-1" complex members:

	deplink(signature.name = "9-1-1", signature = c("RAD9A", "RAD1", "HUS1", "RAD17"))

The results will be output to a local directory (default: root directory) under a folder in name of the designated "signature.name" ("9-1-1" in this case).

Several cutoffs are set by default as below and can be changed by will. Please see the help page for more details (*?deplink*).

	cutoff.freq        = 10
    cutoff.percentile  = 0.2
    cutoff.pvalue      = 0.05
    cutoff.qvalue      = 0.1
    cutoff.diff        = 0.1
    cutoff.fc          = 2

The comparison covers the following features:

* Genomic/epigenetic features
  - [x] Genetic dependency
  - [x] Gene expression
  - [x] Chromatin modification

* Genome instability
  - [x] Genetic mutation
  - [x] COSMIC signature
  - [x] Tumor mutation burden (TMB)
  - [x] Copy number variation (CNV)
  - [x] Microsatellite instability (MSI) 

* Drug sensitivity
  - [x] Drug sensitivity from GDSC data set
  - [x] Drug sensitivity from PRISM data set

* Immune infiltration
  - [x] Immune signature gene (ISG) 

* Stemness
  - [x] mRNA stemness index (mRNAsi)
  - [x] Epithelial–mesenchymal transition (EMT) 

* Misc.
  - [x] Cancer type
  - [x] Hallmark signature


<!--
## :orange_book: What is Programmed Ribosomal Frameshifting (PRF)?

<a href="https://doi.org/10.1016/j.febslet.2013.03.002"><img src="docs/Structural-diversity.png" height="400" align="center" /></a>
- From [*Mauger et al., 2013, FEBS Letters*](https://doi.org/10.1016/j.febslet.2013.03.002)

Ribosomal frameshifting, also known as translational frameshifting or translational recoding, is a biological phenomenon 
that occurs during translation that results in the production of multiple, unique proteins from a single mRNA. 
The process can be programmed by the nucleotide sequence of the mRNA and is sometimes affected by the secondary, 3-dimensional mRNA structure.
It has been described mainly in viruses (especially retroviruses), retrotransposons and bacterial insertion elements, and also in some cellular genes.

For details, please visit [Ribosomal frameshift](https://en.wikipedia.org/wiki/Ribosomal_frameshift).
-->

## :pencil2: Authors

Xiao CHEN, PhD

Herbert Irving Comprehensive Cancer Center, Columbia University Medical Center, New York

<https://www.researchgate.net/profile/Xiao_Chen126>

<!-- [![Twitter](https://img.shields.io/twitter/url/http/shields.io.svg?style=social&logo=twitter)](https://twitter.com/intent/tweet?hashtags=deplink&url=https://github.com/seanchen607/deplink&screen_name=SC607) -->

If you use [deplink](https://github.com/seanchen607/deplink) in
published research, please cite the most appropriate paper(s) from this
list:

1.  **X Chen**, J McGuire, F Zhu, X Xu, Y Li, D Karagiannis, R Dalla-Favera, A Ciccia, J Amengual & C Lu (2020). 
    Harnessing genetic dependency correlation network to reveal chromatin vulnerability in cancer.
    ***In preparation***.
