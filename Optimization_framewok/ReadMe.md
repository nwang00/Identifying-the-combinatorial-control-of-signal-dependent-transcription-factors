# README

## Using optimization framework to identify GRS from perturbation response data
To run the error optimization program, user needs to provide: 1. TF activities, 2. gene expression, 3. data uncertainty of gene expression datasets, 3. file include all the parameters.

### Input: TF activities
Quantified timecourse TF activties read from TF activity file (`./TF_activity_sample.txt` in example). 

Here is the example of TF activties file:
```
NFkB_WT_LipidA	NFkB_TRIF-/-_LipidA	NFkB_MYD88_LipidA	NFkB_JNK-/-/ERK-/-_LipidA	NFkB_Pam3CSK4	IRF3_WT_LipidA	IRF3_TRIF-/-_LipidA	IRF3_MYD88_LipidA	IRF3_JNK-/-/ERK-/-_LipidA	IRF3_Pam3CSK4	MAPK_WT_LipidA	MAPK_TRIF-/-_LipidA	MAPK_MYD88_LipidA	MAPK_JNK-/-/ERK-/-_LipidA	MAPK_Pam3CSK4
0	0	0	0	0	0	0	0	0	0	0.075009871	0.101798537	0.096549387	0	0.084542562
0.127167993	0.60724066	0.122866609	0.127167993	0.772680261	0	0	0.23560308	0	0	1	0.923162118	0.145273963	0	1.84097791
0.747983511	0.965219213	0.33005213	0.747983511	0.800275423	0.14765321	0	1.08732367	0.14765321	0	0.823146158	0.922963564	1.040252829	0	0.651769381
1	0.338405266	0.807611544	1	0.603882419	1	0	1.973889	1	0	0.209820923	0.153563652	0.192755669	0	0.121685025
```
Columns corresponds to different perturbation conditions. Rows corresponds to multiple time points.

### Input: gene expression
Quantified timecourse RNA-seq data is read from gene expression file (`./Gene_exp_fit_sample.txt` in example). Here gene expression is quantified by RPKM. 

Here is the example of gene expression file:
```
WT1_0.bam	WT1_15.bam	WT1_30.bam	WT1_60.bam	TRIF_0_GGCTACA.bam	TRIF_15_TTAGGCA..bam	TRIF_30_GTGAAAC...bam	TRIF_60_GATCAGA..bam	MyD88_0_ATCACGA.bam	MyD88_15_TTAGGCA.bam	MyD88_30_ACTTGAA.bam	MyD88_60_GATCAGA.bam	BP_0_rep_1.bam	BP_15_rep_1.bam	BP_30_rep_1.bam	BP_60_rep_1.bam	X.PAM_0._.GTTTCGG...bam	X.PAM_15._.CGTACGT...bam	X.PAM_30._.GAGTGGA...bam	X.PAM_60._.ACTGATA...bam
0.115608927401906	0.223468891811705	1.38146665245674	63.2257046056402	0	0.0428447960259731	0.0695972425123051	0.828971755992528	0.0945442383439439	0.104637294282298	1.32801209551003	47.0975965477427	0.0592353381262199	0.100700646560097	0.919096620814375	56.0164079827999	0.0107806362342238	0.0549480780313137	0.459844748424772	3.62261737578923
```
The gene expression file include one gene expression pattern under multiple perturbations for optimization. Each number represents mean of multiple replicates as gene expression pattern for optimization.

Note:
Please make sure perturbation conditions and timepoints are consistent between TF activities and gene expression.

### Input: data uncertainty of gene expression datasets
Estimated data uncertainty for timecourse RNA-seq data is read from data uncertainty file (`./Gene_exp_uncertainty.txt` in example). This result can be obtained from error model program. User needs to cpmpose multiple perturbations together into the data uncertainty file. 

Here is the example of data uncertainty file:
```
Uncertainty_class	WT_LipidA	TRIF-/-_LipidA	MYD88_LipidA	JNK-/-/ERK-/-_LipidA	Pam3CSK4
Basal_variance	0.0253	0.0253	0.0328	0.07016	0.0253
Inducible_variance_parameters(x0)	0	0	0	0	0
Inducible_variance_parameters(x1)	0	0	0	0	0
Inducible_variance_parameters(x2)	3.50E-02	3.50E-02	4.18E-02	9.69E-02	3.50E-02
Temporal_varaince	1.56	1.56	3.06	5.63E+01	1.56
```
Columns corresponds to different perturbation conditions. Rows corresponds to different parameters to represent data uncertainty.

Note:
Please make sure perturbation conditions are consistent with previous TF activities and gene expression.


### Input: file include all the parameters
All the parameters can be changed in the file called `Parameters_optimization.txt`. Here, user can specify timecourse for RNA-seq data, numer of replicates of each time points, and file path of gene expression file, through the following parameters: `Time`, `Number_of_replicates`, `Path_of_gene_expression_file` in the parameter file.

Here is the example of parameter file:
```
Path_of_TF_activities_file	./TF_activity_sample.txt
Path_of_gene_exp_file	./Gene_exp_fit_sample.txt
Path_of_gene_exp_uncertainty_file	./Gene_exp_uncertainty.txt
Time_TF	0,15,30,60
Gene_exp_TF	0,15,30,60
Num_of_perturbations 5
```

Notes:
Make sure timecourse is ordered from early to the late time points.

### Command to run the program
To run the program, simply change the directory to error model folder (`./Error_Model`) and run the following command in terminal. Path of the parameter file (`./Parameters_ErrorModel.txt` in example) can be replaced by user's parameter file path.
```
Rscript Error_Model_program.R ./Parameters_ErrorModel.txt
```

### Results
Inferring data uncertainty for timecourse RNA-seq takes time. For the tested example, user should explect to take couple of hours to obtain the results.

After the completion of running the program, the program will return the inferred data uncertainty to a output text file (`Inferred_uncertainly_level.txt`).
```
Basal_variance 12.5754699109369
Inducible_variance_parameters(x0) 2e-20
Inducible_variance_parameters(x1) 0.00235705950636783
Inducible_variance_parameters(x2) 0.269338232494641
Temporal_varaince 1.56378813855195e-11
```

In this file, user can see the inferred data uncertainty represented by multiple parameters including `Basal_variance`, `Inducible_variance_parameters(x0)`, `Inducible_variance_parameters(x1)`, `Inducible_variance_parameters(x2)`, `Temporal_varaince`.
