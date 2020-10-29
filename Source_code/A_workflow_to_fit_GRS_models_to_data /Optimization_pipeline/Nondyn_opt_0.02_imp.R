# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics")
rm(list=ls())
# library(foreach)
# library(iterators)
# library(parallel)
# library(plyr)
# library(doMC)
library("MASS")
library(deSolve)
library(signal)
source("Nondyn_opt_newbasalnew.R")
load("Nondyn_bio_V0.01_B5T5_rep1000_0.02.RData")
load("Prior_TV_bio_V0.01_B5T5_rep1000_0.02_imp.RData")
var_est<-Prior_TV_bio_V0.01_B5T5_rep1000_0.02_imp$var_est_n
var_est_g<-Prior_TV_bio_V0.01_B5T5_rep1000_0.02_imp$var_est_g
Nondyn_V0.01_TB5_TT5_rep1000_p<-cbind(Nondyn_bio_V0.01_B5T5_rep1000_0.02[,1:104],Nondyn_bio_V0.01_B5T5_rep1000_0.02[,105:208])

# doMC::registerDoMC(cores=1)
args<- as.numeric(commandArgs(trailingOnly=TRUE))
# args<-13
run_id<-rownames(Nondyn_V0.01_TB5_TT5_rep1000_p)[args]
name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
#####################################################
#####################################################True TF activities
Basal_level=0
Medium_strength<-0.5
input_list_raw<-as.data.frame(expand.grid(c(1,Medium_strength,Basal_level),c(1,Medium_strength,Basal_level),c(1,Medium_strength,Basal_level)))[-27,]
####change deltaN
deltaN<-1
set.seed(1)
input_error<-matrix(rnorm(3*26,mean=0,sd=deltaN*0.02),nrow=26) 
input_list<-input_list_raw*(1+input_error)
#####################################################
#####################################################
#####################################################Function test
#####################################################
#############Error for RNA
####likelihood of RNA dynamics fitting
slope_ind<-round(2*(var_est[4])^0.5)
if(slope_ind>60) slope_ind<-60
if(slope_ind<1) slope_ind<-1

caRNA_time<-c(0,0,0,(15-slope_ind),15,(15+slope_ind),(30-slope_ind),30,(30+slope_ind),(60-slope_ind),60,(60+slope_ind))
caRNA_time[caRNA_time<0]<-0
caRNA_time[caRNA_time>80]<-80
slope_interval_raw<-c(1,(caRNA_time[6]-caRNA_time[4]),(caRNA_time[9]-caRNA_time[7]),(caRNA_time[12]-caRNA_time[10]))
slope_interval<-rep(rep(slope_interval_raw,each=3),26)
###caRNA_time<-c(0,0,2,14,15,16,29,30,31,58,60,60)
######
# ind_data<-c(2,5,8,11)+rep(seq(0,12*25,12),each=4)
######26 data sets*4 time poitns for each
ind_stdv<-rep(1:104,each=3)
data_ind<-seq(2,104*3,3)
data_ind_full<-rep(seq(2,104*3,3),2)
exp_data<-as.numeric(Nondyn_V0.01_TB5_TT5_rep1000_p[args,])

ind_fun<-function(x){
  return((x<=3))
}
ind_fun_L<-function(x){
  return((x>3))
}

# model_data<-rep(exp_data[1:104],each=3)
#####check
target_fun<-function(model_data){
  f_mean<-var_est[1]+var_est[2]*model_data[data_ind]+var_est[3]*(model_data[data_ind])^2
  slope_list<-tapply(model_data/slope_interval,ind_stdv,function(x) (x[3]-x[1]))
  var_RNA<-rep(((f_mean+(slope_list)^2*var_est[4])*ind_fun_L(model_data[data_ind])+var_est_g),2)
  
  return(sum(-(1-ind_fun(model_data[data_ind_full]))*log(dnorm(exp_data,mean=model_data[data_ind_full],sd=sqrt(var_RNA))+10^(-20))-
               ind_fun(model_data[data_ind_full])*log(dgamma(exp_data, shape=(model_data[data_ind_full])^2/(var_RNA),scale = (var_RNA)/model_data[data_ind_full])+10^(-20))))
}
###################################
###################################
sample_size<-300
K1p<-replicate(8,runif(sample_size, min=-2, max=2))
K2p<-replicate(8,runif(sample_size, min=-2, max=2))
K3p<-replicate(8,runif(sample_size, min=-2, max=2))
K4p<-replicate(8,runif(sample_size, min=-2, max=2))
K5p<-replicate(8,runif(sample_size, min=-2, max=2))
K6p<-replicate(8,runif(sample_size, min=0, max=1))

for(GRN_ind in 1:8){
  tmp_opt<-vector()
  for(i in 1:sample_size) {
    tmp_opt<-c(tmp_opt,optim(c(K1p[i,GRN_ind],K2p[i,GRN_ind],K3p[i,GRN_ind],K4p[i,GRN_ind],K5p[i,GRN_ind],K6p[i,GRN_ind]),method="L-BFGS-B",P1_fun_GRNs,hessian=TRUE))
  } 
  assign(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""),tmp_opt)
  save(list=paste(name_id,"_nondyn_GRN",GRN_ind,sep=""),file = paste(name_id,"_nondyn_Noise_",GRN_ind,".RData",sep=""))
} 
###################################
###################################saving data
list<-vector()
for(args in 1:93){
  run_id<-rownames(Nondyn_V0.01_TB5_TT5_rep1000_p)[as.numeric(args)]
  name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
  for(GRN_ind in 1:8){
    load(paste(name_id,"_nondyn_Noise_",GRN_ind,".RData",sep=""))
    list=c(list,paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))
  }
} 
save(list=list,file = "Nondyn_V0.01_TB5_TT5_0.02_imp.RData")