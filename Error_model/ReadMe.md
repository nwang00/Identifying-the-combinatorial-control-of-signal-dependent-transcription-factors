# README

## Using error model to quantify data uncertainty for timecourse RNA-seq data
To run the error model program, user has to provide both quantified timecourse RNA-seq data and parameter file. 

### Timecourse RNA-seq data as input file
Quantified timecourse RNA-seq data is read from gene expression file ('./Gene_exp_sample.txt' in example). Here gene expression is quantified by RPKM. 

Here is the example of gene expression file:
```
GENE_ID WT1_0.bam WT2_0.bam WT1_15.bam WT2_15.bam WT1_30.bam WT2_30.bam WT1_60.bam WT2_60.bam
NM_001011874 0.000234381819243227 7.69917918186358e-05 0.000275367735930668 7.11419379407099e-05 0 0 0 6.8277908660207e-05
NM_001195662 0.00150788421744375 0.000891581414614788 0.000354313030203234 0.000549225903777989 0.000706863522258849 0.000754300542186464 0.000617314953862357 0.000263557594700512
NM_011283 0.00424863064751238 0 0.00249579042860198 0 0 0.00354220874951248 0 0
```
The gene expression file start from Gene ID, then followed by timecourse RNA-seq data with multiple replicates.


### Alter parameters in the parameter file
All the parameters can be changed in the file called 'Parameters_ErrorModel.txt'. Here, user can specify timecourse for RNA-seq data, numer of replicates of each time points, and file path of gene expression file, through the following parameters: 'Time', 'Number_of_replicates', 'Path_of_gene_expression_file' in the parameter file.

Here is the example of parameter file:
```
Time	0,15,30,60
Number_of_replicates	2
Path_of_gene_expression_file	./Gene_exp_sample.txt
```

Notes:
Make sure timecourse is ordered from early to the late time points.

### Command to run the program
To run the program, simply change the directory to error model folder ("./Error_Model") and run the following command in terminal. Path of the parameter file ("./Parameters_ErrorModel.txt" in example) can be replaced by user's parameter file path.
```
Rscript Error_Model_program.R ./Parameters_ErrorModel.txt
```

### Interpret results

