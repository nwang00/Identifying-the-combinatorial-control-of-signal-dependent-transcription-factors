# README

## Guidelines for designing perturbation studies
To examine identifiability of different GRSs given TF dynamics, user could use Perturbation_analysis_final.R for the analysis. 

### Loading of data
Before running of the downstream analysis, the script will load the gene expression data for each GRS and data uncertainties (`Load Expression Data` section in the script). User could customize the gene expression data for downstream analysis (documented below).

### Alter TF dynamics and gene expression data
To alter TF dynamics and gene expression data for the testing and downstream analysis, user could do so by using script `Nondyn_noise_generate_0.02.R` in the subfolder `Error Model Benchmarking` under `An error model to quantify data uncertainty in value and time`. The details of altering TF dynamics and gene expression are documented in the read me file in `An error model to quantify data uncertainty in value and time`.

### Alter data uncertainty level
To alter the data uncertainty level of gene expression data, user could do so by using the same script mentioned above in `Nondyn_noise_generate_0.02.R` under the folder `Error Model Benchmarking`.

### Downstream analysis
To examine identifiability of GRS from customized gene expression data, user only needs to replace `Nondyn_simulate_0.02.mean` with new gene expression data, and `var_RNA` with new data uncertainty level. The user could automatically generate the figures for downstream analysis by running the rest of the script.
