# optimalCytosolicDensity_linearModel
This repository stores the MATLAB source code of the 'linear model' and also the 'parallel model' in the following publication:

Pang, T. Y., & Lercher, M. J. (2020). Optimal density of bacterial cells. Cold Spring Harbor Laboratory. https://doi.org/10.1101/2020.11.18.388744

See Fig. 2 of the manuscript for an illustration of the model.

The three folders correspond to the 3 versions of the model: (1) main model, (2) diffusion limited case, (3) transition state limited case. See Fig. 1 of the manuscript for an explanation to the different limited cases.

Each folder has the main script files 'script1_simulation.m' and 'script2_make_plots.m'. Script 1 performs the simulation and writes the intermediateData_* files that store the simulation results. Script 2 reads the intermediateData_* files and make the plots (Fig. 2 and Fig. 7 of the manuscript). The other *.m files contains functions called by the main scripts.

All the MATLAB scripts have been tested in MATLAB version 9.10.0.1602886 (R2021a).
