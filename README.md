
<b> About </b>

scDDboost is an R package to analyze changes in the distribution of single-cell expression data between two experimental conditions.
Compared to other methods that assess differential expression, scDDboost benefits uniquely from information conveyed by the clustering of cells into cellular subtypes.  Through a novel empirical Bayesian formulation it calculates gene-specific posterior probabilities that 
the marginal expression distribution is the same (or different) between the two conditions.  The implementation in scDDboost treats gene-level expression data within each condition as a mixture of negative binomial distributions.  


<b> Installation </b>

To install the R package:
```R
# install.packages("devtools")
devtools::install_github("wiscstatman/scDDboost")
```
A tutorial and examples can be found at Rpackage/vignette/ 

<b> Paper </b>

Ma, X., Korthauer, K., Kendziorski, C., and Newton, M. A. (2019). A Compositional Model To Assess Expression Changes From Single-Cell RNA-Seq Data. 
The Annals of Applied Statistics 15, no. 2 (2021): 880-901. https://doi.org/10.1214/20-AOAS1423 .  Earlier version, <a href="https://www.biorxiv.org/content/10.1101/655795v1.abstract"> bioRxiv 655795 </a>

