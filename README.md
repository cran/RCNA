# Robust Copy Number Alteration Detection (RCNA)

RCNA is an R package used to detect the copy number alterations from the targeted panel sequencing data without matched controls.

# Introduction

Copy number alteration (CNA) is a kind of mutation widely observed across a variety of tumor types. 
Recurrent CNAs in cancer cells have been identified as an important factor in cancer development, 
and somatic CNAs promote the development and progression of cancers by altering the gene expression levels.
Sequencing of a panel of targeted genes is common for cancer research given the relatively low cost compared with whole genome sequencing and whole exome sequencing. 
However, detecting CNAs from targeted sequencing remains a challenge due to the limited breadth of the coverage and capture bias in different targeted regions.
Here we describe a new method for robust method for copy number alteration (CNA) detection that is specially designed to detect CNAs from targeted sequencing of tumor samples without matched controls. 
In this method, sequencing coverage bias caused by GC content were corrected before CNA detection to eliminate the different performance on regions from different genes. 
This method was based on sequencing coverage data at base level that enable it to detect CNAs in small scale (e.g. exon level), 
and can be applied to the targeted sequencing dataset with a small number of gene panels. 

# Installation

*For now this package is only offered through the GitHub page, pending further acceptance from CRAN or other repositories*

You may install RCNA from CRAN using the following command:

```
install.packages("RCNA")
```

You may also install `RCNA` using `devtools`:
```
devtools::install_github("repo")`
```
# Tutorial

This package contains tools to 

  - Estimate and correct for GC content bias in sequencing coverage. 
  
  - Calculate normal karyotype range (CI) for CNA estimates
  
  - Estimate copy number alteration (CNA) and report significant alterations.

In order to run RCNA you must have a **sequence file** that specifies position, base, and coverage, as well as an **annotation file** containing genes and their associated genomic intervals. Examples of these files can be found in the vignette.

The general workflow typically looks like 

```
correct_gc_bias() -> estimate_nkr() -> estimate_feature_score()
```

This series of commands will run the RCNA workflow for determining significant CNA after correcting for GC-content bias in coverage. The resulting tables will be in **tab.gz** format and can be found in their respective folders. For more details on the proper execution of the workflow, see each individual functions' help (`?<function-name>`) and the vignette included in this package.
