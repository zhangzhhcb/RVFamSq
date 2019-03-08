# RVFamSq
Rare Variant-Family-Based Score Test for Quantitative Traits  

Description:

The RVFamSq package provides an efficient approach to examine the association between the region-based rare variants and the quantitative traits in family-based data. The RVFamsSq can be broadly applied to diverse pedigrees with members missing sequence data. In addition, qualitative and quantitative covariates, e.g., age, sex, and body mass index, can be flexibly included. The speed of the package is significantly optimized to analyze large-scale data sets.  

Installation:

1.To install the RVFamSq package using the following command:
  install.packages('path/to/RVFamSq_0.1.0.tar.gz', repos = NULL)
  The RVFamSq_0.1.0.tar.gz is the source files that can be downloaded from https://github.com/zhangzhhcb/RVFamSq/tree/master.
2. The depent R-packages used in RVFamSq can be downloaded and installed automatically from CRAN repositories using:
  install.packages(c('bbmle', 'mvtnorm', 'rlist'))
