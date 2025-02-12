# hlacoloc
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13989351.svg)](https://doi.org/10.5281/zenodo.13989351)

Colocalization for the HLA region. See the vignette pdf or instructions below for information on how to make the method work.

You can install this package to your cluster using ***install_github("DrGBL/hlacoloc")*** with either ***devtools*** or ***remotes***.

More information can be found on the preprint available here: https://www.medrxiv.org/content/10.1101/2024.11.05.24316783v1

## How to cite HLA-coloc
To come.

## General steps performed by HLA-coloc

HLA-coloc performs colocalization in two steps, as described in the manuscript. These are:

1. HLA gene causal signature: this uses SuSiE to assign a posterior inclusion probability (PIP) to each HLA allele.

2. HLA PIP Bayesian regression: this uses stanR to perform Bayesian regression of the PIPs in each genes separately.

## Required inputs

The following dataframes are necessary to run it (for each of the two phenotypes):

1. The HLA association summary statistics. This is a dataframe with the following columns:
  - `Name`: the name of the allele. These must be given in IMGT-HLA format, but can be of any field resolution though (as long as it's the same for both phenotypes).
  - `z`: the z-scores of the associations.
  - `beta`: the beta of the associations (if `z` is not provided).
  - `se`: the standard errors (if `z` is not provided).
  - `N`: the sample sizes.

2. The LD matrix (can be either in dataframe or matrix format). This is a square matrix of R coefficient for the HLA alleles above. Importantly, the order of the alleles needs to be the same in the summary statistics dataframe and in the LD matrix here.

## Example

First, load the data:
```{r}
data("ebna","ms","r_ebna","r_ms")
```

