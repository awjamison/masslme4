# masslme4

masslme4 makes it easier to fit a large batch of [lme4](https://github.com/lme4/lme4) models. It supports parallelization, working with large datasets, and enhanced condition handling. This package was originally designed for mass-univariate analysis of neuroimaging and EEG/MEG data where many models are fit independently (e.g. to each voxel or sensor) and a statistic (e.g. a cluster-level statistic) is computed by combining information in the different models. 

## Release notes

This is an alpha version (0.0.0.9000) and all parts are subject to change, including the API.

## Installation

Using the R command line:

```R
install.packages("devtools")
library(devtools)
devtools::install_github("awjamison/masslme4")
```

Or clone the repository:

```
$ git clone https://github.com/awjamison/masslme4
``` 
