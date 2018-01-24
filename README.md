# DECENT
Differential Expression with Capture Efficiency adjustmeNT for single-cell RNA-seq data

C. Ye, TP. Speed, A. Salim (2017) DECENT: Differential Expression with Capture Efficiency AdjustmeNT for Single-Cell RNA-seq Data. bioRxiv. https://doi.org/10.1101/225177

## News
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
de.table <- decent(data.obs = sim$data.obs, cell.type = sim$cell.type, use.spikes = T, 
                   spikes = sim$ercc.obs, spike.conc = sim$ercc.true)

# DECENT without spike-ins
de.table <- decent(data.obs = sim$data.obs, cell.type = sim$cell.type)

# Ground truth can be found in the DE.gene vector.
```
The output object of DE model, no-DE model and LRT will be saved in the working directory (```dir``` argument) as ```decent.DE.rds```, ```decent.noDE.rds``` and ```decent.lrt.rds```. A data frame containing the DE results is returned by the function.

Note that the LRT step also involves optimization of parameters and is currently the bottleneck step.

The function is by default run in parallel using all cores. Specify the number of cores to use by changing the argument ```n.cores```. Use option ```parallel = F``` to run on single core.

By default, cell size factors are estimated using MLE. In some cases, TMM (setting ```normalize = 'TMM'```) gives more accurate estimates.
