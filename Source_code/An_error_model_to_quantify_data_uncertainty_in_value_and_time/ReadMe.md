# README

## Using error model to quantify data uncertainty for timecourse RNA-seq data
To run the error model program, user has to provide both quantified timecourse RNA-seq data and parameter file. 

### Simulate gene expression data with given data uncertainty
All the parameters can be changed in the file called `Parameters_ErrorModel.txt`. Here, user can specify timecourse for RNA-seq data, number of replicates of each time points, and file path of gene expression file, through the following parameters: `Time`, `Number_of_replicates`, `Path_of_gene_expression_file` in the parameter file.

Here is the example of parameter file:
```
Time	0,15,30,60
Number_of_replicates	2
Path_of_gene_expression_file	./Gene_exp_sample.txt
```

Notes:
Make sure timecourse is ordered from early to the late time points.

### Data uncertainty estimation
To run the program, simply change the directory to error model folder (`./Error_Model`) and run the following command in terminal. Path of the parameter file (`./Parameters_ErrorModel.txt` in example) can be replaced by user's parameter file path.
```
Rscript Error_Model_program.R ./Parameters_ErrorModel.txt
```

### Testing of different data uncertainty levels
