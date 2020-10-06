# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Exp_test")
rm(list=ls())
library("MASS")
library(deSolve)
library(signal)

source("Dyn_fit_equa_V3.R")
Gene_ind<- as.numeric(commandArgs(trailingOnly=TRUE))
#####################################################
#####################################################
load("Experimental_data.RData")
int_length<-10^-3
Basal_level<-0

exp_data<-as.numeric(Exp_data_mean[Gene_ind,])
##add sudo counts for zero counts
exp_data[exp_data==0]<-0.1

exp_data_WT_max<-max(exp_data[1:4])
exp_data_basal<-rep(exp_data[seq(1,20,4)],each=64)
num_cond<-length(exp_data)/4
#####################################################
load("weak_threshold.RData")
####################################################
#############Error for RNA
####################################################
load("WT_LipidA_GeneExp.Prior.RData")
load("MAPKi_LipidA_GeneExp.Prior.RData")
load("MYD88_LipidA_GeneExp.Prior.RData")

var_est_g_list<-rbind(WT_LipidA_GeneExp.Prior$var_est_g,WT_LipidA_GeneExp.Prior$var_est_g,MYD88_LipidA_GeneExp.Prior$var_est_g,MAPKi_LipidA_GeneExp.Prior$var_est_g,WT_LipidA_GeneExp.Prior$var_est_g)
var_est_n_list<-rbind(WT_LipidA_GeneExp.Prior$var_est_n,WT_LipidA_GeneExp.Prior$var_est_n,MYD88_LipidA_GeneExp.Prior$var_est_n,MAPKi_LipidA_GeneExp.Prior$var_est_n,WT_LipidA_GeneExp.Prior$var_est_n)
###############################calculate slope#################
slope_ind<-round(2*(c(0,var_est_n_list[,4]))^0.5)
slope_ind[slope_ind>60]<-60
slope_ind[slope_ind<1]<-1

caRNA_time_org<-c(c(0,15,30,60),as.vector(t(rbind(c(15,15,30,30,60,60),c(15,15,30,30,60,60))[rep(1,num_cond*2),]+c(-slope_ind,slope_ind))))
caRNA_time_org[caRNA_time_org<0]<-0
caRNA_time_org[caRNA_time_org>80]<-80
slope_interval_g<-(caRNA_time_org[seq(6,length(caRNA_time_org),2)]-caRNA_time_org[seq(5,length(caRNA_time_org),2)])[1:(num_cond*3)]
slope_interval_n<-(caRNA_time_org[seq(6,length(caRNA_time_org),2)]-caRNA_time_org[seq(5,length(caRNA_time_org),2)])[(num_cond*3+1):(num_cond*6)]
caRNA_time<-sort(caRNA_time_org)
back_order<-rep(order(order(caRNA_time_org)),num_cond)+rep(seq(0,64*(num_cond-1),64),each=64)
################################################
################################################
data_ind<-rep(c(1:4),num_cond)+rep(seq(0,64*(num_cond-1),64),each=4)
slope_list_n<-rep(0,4*num_cond)
slope_list_g<-rep(0,4*num_cond)

slop_list_ind<-c(1:(4*num_cond))[-seq(1,4*num_cond,4)]

slop_cal_ind1_g<-c(5:10,11:16+64,17:22+64*2,23:28+64*3,29:34+64*4)[seq(1,30,2)]
slop_cal_ind2_g<-c(5:10,11:16+64,17:22+64*2,23:28+64*3,29:34+64*4)[seq(2,30,2)]

slop_cal_ind1_n<-(c(5:10,11:16+64,17:22+64*2,23:28+64*3,29:34+64*4)+30)[seq(1,30,2)]
slop_cal_ind2_n<-(c(5:10,11:16+64,17:22+64*2,23:28+64*3,29:34+64*4)+30)[seq(2,30,2)]

var_est_g_cal<-var_est_g_list[rep(c(1:num_cond),each=4),]
var_est_n_cal<-var_est_n_list[rep(c(1:num_cond),each=4),]

ind_fun<-function(x){
  return((x<=3))
}
##############################
target_fun<-function(model_data){
  
  model_data<-model_data+exp_data_basal
  model_data_reorder<-model_data[back_order]
  
  f_mean_n<-var_est_n_cal[,1]+var_est_n_cal[,2]*model_data_reorder[data_ind]+var_est_n_cal[,3]*(model_data_reorder[data_ind])^2
  slope_list_n[slop_list_ind]<-(model_data_reorder[slop_cal_ind2_n]-model_data_reorder[slop_cal_ind1_n])/slope_interval_n
  var_RNA_n<-f_mean_n+(slope_list_n)^2*var_est_n_cal[,4]
  
  var_RNA_g<-var_est_g_cal
  
  return(sum(-(1-ind_fun(model_data_reorder[data_ind]))*log(dnorm(exp_data,mean=model_data_reorder[data_ind],sd=sqrt(var_RNA_n))+10^(-20))-
               ind_fun(model_data_reorder[data_ind])*log(dgamma(exp_data, shape=(model_data_reorder[data_ind])^2/(var_RNA_g),scale = (var_RNA_g)/model_data_reorder[data_ind])+10^(-20))))
}
###################################
set.seed(18)
sample_size_full<-5000
sample_size<-30

K1p_raw<-replicate(8,runif(sample_size_full, min=-2, max=5))
K2p_raw<-replicate(8,runif(sample_size_full, min=-2, max=5))
K3p_raw<-replicate(8,runif(sample_size_full, min=-2, max=5))
K5p<-replicate(8,runif(sample_size, min=-2, max=1))
#########Remove parameters combinations for weak gate output####
for(GRN_ind in 1:8){
  for(i in 1:sample_size_full){
    Gate_output<-get(paste("G_GRN",GRN_ind,sep=""))(c(Kd1=10^(K1p_raw[i,GRN_ind]),Kd2=10^(K2p_raw[i,GRN_ind]),Kd3=10^(K3p_raw[i,GRN_ind]),h1=1,h2=1,h3=1))
    if(Gate_output < weak_threshold) {
      K1p_raw[i,GRN_ind]<-NA
      K2p_raw[i,GRN_ind]<-NA
      K3p_raw[i,GRN_ind]<-NA
    }
  }
}

K1p<-matrix(NA,ncol = 8,nrow=sample_size)
K2p<-matrix(NA,ncol = 8,nrow=sample_size)
K3p<-matrix(NA,ncol = 8,nrow=sample_size)

for(GRN_ind in 1:8){
  K1p[,GRN_ind]<-K1p_raw[which(!is.na(K1p_raw[,GRN_ind]))[1:sample_size],GRN_ind]
  K2p[,GRN_ind]<-K2p_raw[which(!is.na(K2p_raw[,GRN_ind]))[1:sample_size],GRN_ind]
  K3p[,GRN_ind]<-K3p_raw[which(!is.na(K3p_raw[,GRN_ind]))[1:sample_size],GRN_ind]
}
########################################################################
########################################################################
ksyn_cal<-function(logKd1,logKd2,logKd3,logkdeg,logksyn,Basal_level,exp_data_WT_max){
  R_basal<-10^(logksyn)*(get(paste("G_basal",GRN_ind,sep=""))(10^(logKd1),10^(logKd2),10^(logKd3),Basal_level))/10^(logkdeg)
  return(exp_data_WT_max/max(ode(y=c(RNA1=R_basal,RNA2=R_basal,RNA3=R_basal,RNA4=R_basal,RNA5=R_basal),
                                 times = c(0,15,30,60),func =get(paste("P1_GRN",GRN_ind,"_fit",sep="")),parms = c(Kd1=10^(logKd1),Kd2=10^(logKd2),Kd3=10^(logKd3),ksyn=10^(logksyn),kdeg=10^(logkdeg)),hmax=1)[,2]))
}

for(GRN_ind in 1:8){
  tmp_opt<-vector()
  for(i in 1:sample_size) {
    # start_time <- Sys.time()
    Ksyn_sample<-ksyn_cal(logKd1 = K1p[i,GRN_ind],logKd2 = K2p[i,GRN_ind],logKd3 = K3p[i,GRN_ind],logksyn=0,logkdeg=K5p[i],Basal_level = Basal_level,exp_data_WT_max=exp_data_WT_max)
    tmp_opt<-c(tmp_opt,optim(c(K1p[i,GRN_ind],K2p[i,GRN_ind],K3p[i,GRN_ind],log10(Ksyn_sample),K5p[i]),method="L-BFGS-B",lower = c(rep(-2,3),-5,-2),upper = c(rep(5,3),6,1),P1_fun_GRNs,hessian=F))
    # print(c(i,(Sys.time()-start_time),tmp_opt[i*5-3]$value,tmp_opt[i*5-4]$par))
  } 
  assign(paste("dyn_Gene",Gene_ind,"_Logic",GRN_ind,sep=""),tmp_opt)
  save(list=paste("dyn_Gene",Gene_ind,"_Logic",GRN_ind,sep=""),file = paste("dyn_Gene",Gene_ind,"_Logic",GRN_ind,".RData",sep=""))
}


##############################################################################
####################save files#######################################
##############################################################################
list<-vector()
for(GRN_ind in 1:8){
  load(paste("dyn_Gene",Gene_ind,"_Logic",GRN_ind,".RData",sep=""))
  list=c(list,paste("dyn_Gene",Gene_ind,"_Logic",GRN_ind,sep=""))
}
save(list=list,file = paste("dyn_Gene",Gene_ind,".RData",sep=""))