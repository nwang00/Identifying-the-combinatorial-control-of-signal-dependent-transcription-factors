# README

## Guidelines for designing perturbation studies
To examine identifiability of different GRSs given TF dynamics, user could use Perturbation_analysis_final.R for the analysis. 

### Loading of data
Before running of the analysis, user needs to load the gene expression data for each GRS and data uncertainties (Load Expression Data section in the script).

### Alter TF dynamics and gene expression data
All the parameters can be changed in the file called `Parameters_ErrorModel.txt`. Here, user can specify timecourse for RNA-seq data, number of replicates of each time points, and file path of gene expression file, through the following parameters: `Time`, `Number_of_replicates`, `Path_of_gene_expression_file` in the parameter file.

Here is the example of parameter file:
```
Time	0,15,30,60
Number_of_replicates	2
Path_of_gene_expression_file	./Gene_exp_sample.txt
```

Notes:
Make sure timecourse is ordered from early to the late time points.

### Alter data uncertainty level
To run the program, simply change the directory to error model folder (`./Error_Model`) and run the following command in terminal. Path of the parameter file (`./Parameters_ErrorModel.txt` in example) can be replaced by user's parameter file path.
```
Rscript Error_Model_program.R ./Parameters_ErrorModel.txt
```

### Results
Inferring data uncertainty for timecourse RNA-seq takes time. For the tested example, user should expect to take couple of hours to obtain the results.

After the completion of running the program, the program will return the inferred data uncertainty to a output text file (`Inferred_uncertainly_level.txt`).
```
Basal_variance 12.5754699109369
Inducible_variance_parameters(x0) 2e-20
Inducible_variance_parameters(x1) 0.00235705950636783
Inducible_variance_parameters(x2) 0.269338232494641
Temporal_varaince 1.56378813855195e-11
```

In this file, user can see the inferred data uncertainty represented by multiple parameters including `Basal_variance`, `Inducible_variance_parameters(x0)`, `Inducible_variance_parameters(x1)`, `Inducible_variance_parameters(x2)`, `Temporal_varaince`.
