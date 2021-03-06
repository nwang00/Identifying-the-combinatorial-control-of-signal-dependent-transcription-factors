# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Revise_testing")
rm(list=ls())
library("MASS")
library(deSolve)
library(signal)

source("Dyn_fit_equa_V3.R")

##Ppp1r15a
# Gene_ind<- 16
# var_ind<-1
##Srgn
Gene_ind<- 7
var_ind<-2
#####################################################
#########################Raw variance calculation############################
load("RPKM_Smale.RData")
WT_lipidA_mean_matirx<-(RPKM_Smale[,29:32]+RPKM_Smale[,34:37])/2
Control_MAPKi_lipidA_mean_matirx<-(RPKM_Smale[,c(71,75,77,81)]+RPKM_Smale[,c(70,74,76,80)])/2
norm_factor_raw<-WT_lipidA_mean_matirx/(Control_MAPKi_lipidA_mean_matirx+0.1)
norm_factor<-norm_factor_raw
norm_factor[norm_factor==0]<-1
norm_factor_final<-norm_factor[,c(1,1,2,2,3,3,4,4)]

MAPKi_LipidA_GeneExp_raw<-RPKM_Smale[,c(40,39,44,43,46,45,50,49)]*norm_factor_final

GeneExp_data<-cbind(RPKM_Smale[,c(29,34,30,35,31,36,32,37,100,101,104:109)],MAPKi_LipidA_GeneExp_raw)
###########
Gene_select_id<-which(RPKM_Smale[,2]==c("Ppp1r15a","Srgn"))
Var_Ppp1r15a<-apply(matrix(GeneExp_data[Gene_select_id[1],],nrow=2),2,function(x) var(unlist(x)))
Var_Srgn<-apply(matrix(GeneExp_data[Gene_select_id[2],],nrow=2),2,function(x) var(unlist(x)))
Var_raw<-rbind(Var_Ppp1r15a,Var_Srgn)[var_ind,c(1:4,1:12,1:4)]
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
caRNA_time<-c(0,15,30,60)
data_ind<-rep(c(1:4),num_cond)+rep(seq(0,64*(num_cond-1),64),each=4)

ind_fun<-function(x){
  return((x<=3))
}
##############################
target_fun<-function(model_data){
  model_data[which(model_data==0)]<-0.1

  return(sum(-(1-ind_fun(exp_data))*log(dnorm(model_data,mean=exp_data,sd=sqrt(Var_raw))+10^(-20))-
               ind_fun(exp_data)*log(dgamma(model_data, shape=(exp_data)^2/(Var_raw),scale = (Var_raw)/exp_data)+10^(-20))))
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