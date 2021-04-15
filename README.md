# utoaztecan-analyses

The analyses for Greenhill, Haynie et al. subm. 

Please do not use or cite without permission until this paper is published. 

Contact greenhill@shh.mpg.de for more information


## Contents:

### ./phylogeny/

* ua-covarion-relaxed.xml - BEAST2 XML analysis file.
* ua-covarion-relaxed.log.gz - BEAST2 MCMC log file (gzipped).
* ua-covarion-relaxed.mcct.trees - Maximum Clade Credibility summary tree.
* ua-covarion-relaxed.trees.gz - Posterior probability distribution of trees (gzipped).


### ./ancestral_state/

This folder contains the R code and data necessary for running the ancestral state analysis. 
The [Makefile](https://en.wikipedia.org/wiki/Make_(software)) should work to run all the files
if you run this command from the command line:

> make all

#### Data files:

* EA005.tsv - Recoded data from D-PLACE variable EA005.
* EA028.tsv - Recoded data from D-PLACE variable EA028.
* EA029.tsv - Recoded data from D-PLACE variable EA029.
* EA042.tsv - Recoded data from D-PLACE variable EA042.

#### Code/Analysis files:

* Makefile - `make` file to run all analyses
* get_rates.R - R script for running rates analysis
* plot.R - R script for plotting result
* prettytables.R - R script for tabulating the results
* probabilities.R - R script for summarising the probabilities
* table.R - R script for tabulating the results

#### Result files:

* EA005.result.csv - CSV file containing ancestral state reconstruction estimates.
* EA028.result.csv - CSV file containing ancestral state reconstruction estimates.
* EA029.result.csv - CSV file containing ancestral state reconstruction estimates.
* EA042.result.csv - CSV file containing ancestral state reconstruction estimates.


