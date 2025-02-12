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

```{r}
head(ebna)
#>         Name         beta       se    N
#> 1 G*01:01:01  0.014137500 0.015906 7247
#> 2 G*01:04:04  0.076847400 0.070867 7247
#> 3 G*01:01:02 -0.017568100 0.018972 7247
#> 4 G*01:01:03  0.000931893 0.033494 7247
#> 5   G*01:05N -0.293138000 0.078049 7247
#> 6 G*01:06:01  0.002110380 0.035066 7247
r_ebna[1:5,1:5]
#>            G*01:01:01 G*01:04:04 G*01:01:02 G*01:01:03   G*01:05N
#> G*01:01:01   1.000000 -0.1573650 -0.5404320 -0.2478420 -0.0925080
#> G*01:04:04  -0.157365  1.0000000 -0.0359779 -0.0317737 -0.0171784
#> G*01:01:02  -0.540432 -0.0359779  1.0000000 -0.1419840 -0.0672568
#> G*01:01:03  -0.247842 -0.0317737 -0.1419840  1.0000000 -0.0257557
#> G*01:05N    -0.092508 -0.0171784 -0.0672568 -0.0257557  1.0000000
```
