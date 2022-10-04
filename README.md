# optimalCytosolicDensity_linearModel

## A Model to Investigate the Optimal Cytosolic Density of a Cell

This repository stores the MATLAB source code of the `linear model` (and also the `parallel model`), which investigates how a bacterial cell optimizes its cytosolic density.

## Full Documentation of the Model

Pang, T. Y., & Lercher, M. J. (2020). Optimal density of bacterial cells. Cold Spring Harbor Laboratory. https://doi.org/10.1101/2020.11.18.388744

See Fig. 2 of the manuscript for an illustration of the model.

## Description of the Files

The folders correspond to different versions of the model: (1) main model, (2) diffusion limited case, (3) transition state limited case. See Fig. 1 of the manuscript for an explanation to the different limited cases. 

Each folder has the main script files `script1_simulation.m` and `script2_make_plots.m`. Script 1 performs the simulation and writes the `intermediateData_*` files that store the simulation results. Script 2 reads the `intermediateData_*` files and make the plots (Fig. 3 and Fig. 7 of the manuscript). The other `*.m` files contains functions called by the main scripts.

In addition, the main model is modified to account for crowding using Vasquez's formulation (model 1a); this model is plotted in Fig. S1. The parameter θ within the main model is also modified in a way that the smaller catalysts and metabolites are more diffusive than the larger ones (model 1b); this model is plotted in Fig. 3.

All the MATLAB scripts have been tested in MATLAB version 9.10.0.1602886 (R2021a).
