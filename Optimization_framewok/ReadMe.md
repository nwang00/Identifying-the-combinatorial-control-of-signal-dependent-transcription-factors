# README

## Using optimization framework to identify GRS from perturbation response data
To run the optimization program, user needs to provide: 1. TF activities, 2. gene expression, 3. data uncertainty of gene expression datasets, 4. parameter file that include all the parameters.

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
The gene expression file include gene expression pattern for a single gene to optimize. Each number represents mean of multiple replicates from quantified gene expression file for optimization.

Note:
Please make sure perturbation conditions and timepoints are consistent between TF activities and gene expression.

### Input: data uncertainty of gene expression datasets
Estimated data uncertainty for timecourse RNA-seq data is read from data uncertainty file (`./Gene_exp_uncertainty.txt` in example). This result can be obtained from error model program. User needs to compose multiple perturbations together into the data uncertainty file. 

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
All the parameters can be changed in the file called `Parameters_optimization.txt`. Here, user can specify file paths to TF activities, gene expression, data uncertainty of gene expression, timecourse for TF and gene expression, and number of perturbations, through the following parameters: `Path_of_TF_activities_file`, `Path_of_gene_exp_file`, `Path_of_gene_exp_uncertainty_file`, `Time_TF`, `Time_gene_exp`, `Num_of_perturbations` in the parameter file.

Here is the example of parameter file:
```
Path_of_TF_activities_file	./TF_activity_sample.txt
Path_of_gene_exp_file	./Gene_exp_fit_sample.txt
Path_of_gene_exp_uncertainty_file	./Gene_exp_uncertainty.txt
Time_TF	0,15,30,60
Time_gene_exp	0,15,30,60
Num_of_perturbations 5
```

Notes:
Make sure timecourse is ordered from early to the late time points.

### Command to run the program
To run the program, simply change the directory to optimization program folder (`./Optimization_framework`) and run the following command in terminal. Path of the parameter file (`./Parameters_optimization.txt` in example) can be replaced by user's parameter file path.

```
Rscript Optimization_framework.R ./Parameters_optimization.txt
```

### Results
Identify GRS for gene takes time. For the tested example, user should explect to take about an hour to obtain the results.

After the completion of running the program, the program will return identified GRS to a output text file (`Optimization_output.txt`).
```
Negative_Log_Likelihood Kd1(TF1) Kd2(TF2) Kd3(TF3) ksyn kdeg
TF1 AND TF2 AND TF3 46.2971530619891 0.60439190677249 2.7409160279529 2725.71103911315 118.188757607332 0.340846608150139
TF1 OR TF2 OR TF3 49.3883701553363 1e+05 4.58869147142393 4368.40422298028 32.3425646850761 0.0847699781545087
TF1 OR TF2 AND TF3 49.5528542954707 0.207762777601414 4.95675177437921 3845.63886115806 33.35521921855 0.0796416839879395
TF1 AND (TF2 OR TF3) 51.256476425128 0.0634755995923682 1.43323110134463 119.250678047074 7.2230132828024 0.01
TF2 OR TF1 AND TF3 105.492514130795 0.017261825781774 0.653780687433532 2.1616686572685 275.320858293789 0.343573186444975
TF2 AND (TF1 OR TF3) 140.032322631413 1.78129012763029 0.0152855365037 6.09545093216132 67.7313607473608 0.465450595599247
TF3 OR TF1 AND TF2 143.330620158025 1.84515101976985 0.0216450236093997 52.7515460271333 12.8638566652014 0.0755512885401792
TF3 AND (TF1 OR TF2) 151.446811754697 11.1133441356583 0.729256348468756 0.0151527439980954 14.3524253390713 0.170545972120507
```

In this file, user can see the inferred GRS from 8 logic gates (ordered from good fit to poor fit). The output also includes fitness (represented by negative log likelihood), and all the other estimated parameters (`Kd1(TF1)`, `Kd2(TF2)`, `Kd3(TF3)`, `ksyn`, `kdeg`).
