# DECENT
Differential Expression with Capture Efficiency adjustmeNT for single-cell RNA-seq data

C. Ye, TP. Speed, A. Salim (2018) DECENT: Differential Expression with Capture Efficiency AdjustmeNT for Single-Cell RNA-seq Data. bioRxiv. https://doi.org/10.1101/225177

## News
Jun 21, 2019
* Version 1.1.0 released
* Improved starting values for EM algorithm.
* Reduced memory requirement for single imputation function.
* Other minor changes.

Jun 5, 2019
* Version 1.0.0 released
* Improved bound for GQ integration in the E-step and LRT.
* Improved starting values for EM algorithm.
* Improved global tau estimation.
* Other minor changes.


Feb 16, 2019
* Version 0.99.2 released
* Modified cell-specific tau estimation.
* Other minor changes

Aug 15, 2018
* Version 0.99.1 released
* Rho now depends on the mean by a logistic linear model.
* Gaussian quadrature approximation for optimization in LRT.
* Added single imputation function.
* GLM framework for complex design.
* Other corresponding changes.

Feb 6, 2018
* Version 0.2.0 released.
* Rewrite LRT in Rcpp.

Jan 30, 2018
* Version 0.1.2 released.
* Changed LRT starting values.
* Other minor changes.

Jan 15, 2018
* Version 0.1.1 released.
* Imputed data matrix can now be obtained by calling function `getImputed` or set `imputed = TRUE` when calling `decent`.

## Installation

You can install DECENT from github with:

```R
require(devtools)
devtools::install_github("cz-ye/DECENT")
```

## Quick start

Here we use a simulated dataset for demonstration


```R
data("sim")


# DECENT with spike-ins
de.table <- decent(data.obs = sim$data.obs, # UMI count matrix after quality control
                                            # at least > 3% non-zero counts for each cell and > 5 non-zero counts for each gene
                   X = ~as.factor(sim$cell.type), # cell type/group indicator
                   use.spikes = T, 
                   spikes = sim$sp.obs, # observed UMI count
                   spike.conc = sim$sp.true, # nominal molecule count
                   s.imputed = T, # get single imputation expression matrix
                   E.imputed = T, # get mean imputation expression matrix
                   dir = './' # directory to save the fitted models and imputed data matrices.
                   )

# DECENT without spike-ins
de.table <- decent(data.obs = sim$data.obs,
                   X = ~as.factor(sim$cell.type), 
                   use.spikes = F,
                   CE.range = c(0.02, 0.1) # specify the range of the ranked random capture efficiency
                   )

# DECENT with batch dummy variable
batch <- rep(1, length(sim$cell.type))
set.seed(0)
batch[sample.int(length(sim$cell.type), length(sim$cell.type)/2)] <- 2 # randomly split into 2 batches just for demonstration
de.table <- decent(data.obs = sim$data.obs, 
                   X = ~as.factor(sim$cell.type), 
                   W = ~as.factor(batch),
                   use.spikes = T,
                   spikes = sim$sp.obs, spike.conc = sim$sp.true)
                   
# Ground truth can be found in the DE.gene vector.
```
The output object of DE model, no-DE model and LRT will be saved in the working directory (```dir``` argument) as ```decent.DE.rds```, ```decent.noDE.rds``` and ```decent.lrt.rds```. A data frame containing the DE results is returned by the function.

Note that the LRT step also involves optimization of parameters and is currently the bottleneck step.

The function is by default run in parallel using all cores. Specify the number of cores to use by changing the argument ```n.cores```. Use option ```parallel = F``` to run on single core.

By default, cell size factors are estimated using MLE. In some cases, TMM (setting ```normalize = 'TMM'```) gives more accurate estimates.
