# README

### Simulate gene expression data with given data uncertainty
User could simulate customized gene expression data and data uncertainty level by using script `Nondyn_noise_generate_0.02.R`.

#### Alter gene expression data 
In `Nondyn_noise_generate_0.02.R`, user could alter TF dynamics by alter `input_raw`. User could also alter GRS model by changing `P1_GRN` in function script `Nondyn_opt_newbasalnew.R`. 

#### Alter data uncertainty level 
In `Nondyn_noise_generate_0.02.R`, user could alter data uncertainty level for simulated gene expression data by changing error model parameters. The error model parameters are documented below:

GRS model data uncertainty: `para_noise` for Kd, ksyn and kdeg. `k0_noise` for k0

TF activity data uncertainty: `input_error`

Temporal data uncertainty: `Nondyn_active_temp_B5T5_parV_rep1000_N_0.02`

Sequencing data uncertainty: `Nondyn_bio_V0.01_B5T5_rep1000_0.02`

### Data uncertainty estimation
To estimate data uncertainty of simulated data, user can simply run the script `Error_Model_Infer_0.02.R`. As the data uncertainty estimation has been standardized in the folder `Error Model`, we would suggest user simply use the script `Error_Model_program.R` in the folder `Error Model` for data uncertainty estimation.

### Testing of different data uncertainty levels
In folder `Different_error_levels`, we tested error model performance with different data uncertainty levels by altering error model parameters in `Nondyn_noise_generate_0.02.R`. User can examine error model performance by running the scripts in each subfolder.
