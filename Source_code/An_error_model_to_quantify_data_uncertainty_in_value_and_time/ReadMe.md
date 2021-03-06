# README

### Simulate gene expression data with given data uncertainty
User could simulate customized gene expression data and data uncertainty level by using script `Nondyn_noise_generate_0.02.R`.

#### Alter gene expression data 
In `Nondyn_noise_generate_0.02.R`, user could alter TF dynamics by alter `input_raw`. User could also alter GRS model by changing `P1_GRN` in function script `Nondyn_opt_newbasalnew.R`. 

#### Alter data uncertainty level 
In `Nondyn_noise_generate_0.02.R`, user could alter data uncertainty level for simulated gene expression data by changing error model paramters. The error model parameters are documented below:


### Data uncertainty estimation
To run the program, simply change the directory to error model folder (`./Error_Model`) and run the following command in terminal. Path of the parameter file (`./Parameters_ErrorModel.txt` in example) can be replaced by user's parameter file path.
```
Rscript Error_Model_program.R ./Parameters_ErrorModel.txt
```

### Testing of different data uncertainty levels
