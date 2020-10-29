# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Source_code/Guidelines_for_designing_perturbation_studies")
rm(list=ls())
library("MASS")
library(deSolve)
library(signal)
library("ggplot2")
library(dendextend)
#####################Load Expression Data################################
load("Nondyn_active_GT.RData")
load("Nondyn_bio_V0.01_B5T5_rep1000_0.02.RData")

Nondyn_simulate_0.02<-cbind(Nondyn_bio_V0.01_B5T5_rep1000_0.02[,1:104],Nondyn_bio_V0.01_B5T5_rep1000_0.02[,105:208])
Nondyn_simulate_0.02.mean<-(Nondyn_simulate_0.02[,1:104]+Nondyn_simulate_0.02[,105:208])/2
#####################################################
#####################################################
##################################Error Model##################################
#####################################################
#####################################################
ind_fun<-function(x){
  return((x<=3))
}

target_fun<-function(data,exp_ind,model_ind,perturb_ind,var_GT){
  var_RNA<-rep(var_GT[model_ind,perturb_ind],2)
  data_model<-data[model_ind,c(perturb_ind,perturb_ind+104)]
  data_exp<-data[exp_ind,c(perturb_ind,perturb_ind+104)]
  
  return(sum(-(1-ind_fun(data_model))*log(dnorm(data_exp,mean=data_model,sd=sqrt(var_RNA))+10^(-20))-
               ind_fun(data_model)*log(dgamma(data_exp, shape=(data_model)^2/(var_RNA),scale = (var_RNA)/data_model)+10^(-20))))
}
#############################Error estimate########################
load("var_empirical_bio_V0.01_B5T5_rep1000_0.02.RData")
################################################################################################
################################################################################################
perturb_analysis<-function(Nondyn_simulate=Nondyn_simulate,var_empirical_bio=var_empirical_bio,perturb_raw=perturb_raw){
  perturb_ind<-rep(perturb_raw*4,each=4)-rep(seq(3,0,-1),length(perturb_raw))
  dis<-matrix(0,ncol=93,nrow=93)
  for(i in 1:93){
    for(j in 1:93){
      dis[i,j]<-target_fun(data=Nondyn_simulate,perturb_ind=perturb_ind,exp_ind=i,model_ind=j,var_GT=var_empirical_bio)
    }
  }
  iden<-vector()
  for(GRS_ind in 1:93){
    iden<-c(iden,min(dis[GRS_ind,-GRS_ind])-dis[GRS_ind,GRS_ind])
  }
  return(iden)
}

perturb_all<-perturb_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1:26))
perturb_HHL<-perturb_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,3,7,19))
perturb_HLL<-perturb_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,9,21,25))
perturb_HMM<-perturb_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,13,11,5))
perturb_HHM<-perturb_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,2,4,10))
perturb_HHL_HLL<-perturb_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,3,7,19,9,21,25))
perturb_HMM_HHM<-perturb_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,13,11,5,2,4,10))
perturb_HHL_HLL_HMM_HHM<-perturb_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,3,7,19,9,21,25,13,11,5,2,4,10))
########################################################################
#########################################Plot###############################
########################################################################
Iden_plot<-data.frame(Perturb=rep(c("HLL","HHL","HMM","HHM","HLL_HHL","HMM_HHM","HLL_HHL_HMM_HHM","All"),each=93),Iden=c(perturb_HLL,perturb_HHL,perturb_HMM,perturb_HHM,perturb_HHL_HLL,perturb_HMM_HHM,perturb_HHL_HLL_HMM_HHM,perturb_all))
Iden_plot$Perturb<-factor(Iden_plot$Perturb, levels = c("HLL","HHL","HMM","HHM","HLL_HHL","HMM_HHM","HLL_HHL_HMM_HHM","All"))

pp <- ggplot(Iden_plot, aes(x=Perturb,y=Iden)) + geom_boxplot(mapping=aes(fill=Perturb),outlier.shape = NA,color="grey20",alpha = 0.7) 
pp+labs(x = "",y = "")+
  geom_hline(yintercept=0,color = adjustcolor( "orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
  geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),color = "dodgerblue4",alpha = 0.4)+
  scale_fill_manual(values=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(0),face="bold"),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Iden_perturb_analsis_final.pdf",width=7.5,height=4)
########################################################################
############################Replicates##################################
########################################################################
rep_num<-13
Nondyn_simulate_0.02<-cbind(Nondyn_bio_V0.01_B5T5_rep1000_0.02[,1:(104*rep_num)])
########################################################################
#########################Function###############################################
########################################################################
ind_fun<-function(x){
  return((x<=3))
}

target_fun_rep<-function(data,exp_ind,model_ind,perturb_rep,perturb_proc,var_GT){
  var_RNA<-var_GT[model_ind,perturb_rep]
  data_model<-data[model_ind,perturb_proc]
  data_exp<-data[exp_ind,perturb_proc]
  
  return(sum(-(1-ind_fun(data_model))*log(dnorm(data_exp,mean=data_model,sd=sqrt(var_RNA))+10^(-20))-
               ind_fun(data_model)*log(dgamma(data_exp, shape=(data_model)^2/(var_RNA),scale = (var_RNA)/data_model)+10^(-20))))
}
########################################################################
########################################################################
#######################################################################
replicates_analysis<-function(Nondyn_simulate=Nondyn_simulate,var_empirical_bio=var_empirical_bio,perturb_raw=perturb_raw,rep_num){
  # perturb_raw<-c(1,3,7,19)
  perturb_ind<-rep(perturb_raw*4,each=4)-rep(seq(3,0,-1),length(perturb_raw))
  
  perturb_rep<-rep(perturb_ind,rep_num)
  perturb_proc<-perturb_ind+rep(seq(0,104*(rep_num-1),104),each=length(perturb_ind))
  
  dis<-matrix(0,ncol=93,nrow=93)
  for(i in 1:93){
    for(j in 1:93){
      dis[i,j]<-target_fun_rep(data=Nondyn_simulate,exp_ind=i,model_ind=j,perturb_rep=perturb_rep,perturb_proc=perturb_proc,var_GT=var_empirical_bio)
    }
  }
  iden<-vector()
  for(GRS_ind in 1:93){
    iden<-c(iden,min(dis[GRS_ind,-GRS_ind])-dis[GRS_ind,GRS_ind])
  }
  return(iden)
}

replicates_p4<-replicates_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,3,7,19),13)
replicates_p13<-replicates_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,3,7,19,9,21,25,13,11,5,2,4,10),4)
########################################################################
#########################################Plot###############################
########################################################################
Iden_plot<-data.frame(Perturb=rep(c("Rep4","Rep13","All"),each=93),Iden=c(replicates_p4,replicates_p13,perturb_all))
Iden_plot$Perturb<-factor(Iden_plot$Perturb, levels = c("Rep4","Rep13","All"))

pp <- ggplot(Iden_plot, aes(x=Perturb,y=Iden)) + geom_boxplot(mapping=aes(fill=Perturb),outlier.shape = NA,color="grey20",alpha = 0.7) 
pp+labs(x = "",y = "")+
  geom_hline(yintercept=0,color = adjustcolor( "orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
  geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),color = "dodgerblue4",alpha = 0.4)+
  scale_fill_manual(values=c("dodgerblue4", "dodgerblue4", "dodgerblue4"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(0),face="bold"),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Iden_replicates_analsis_final_revise.pdf",width=5,height=4)
########################################################################
############################Replicates 8##################################
########################################################################
replicates_s1<-replicates_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1),8)
replicates_s2<-replicates_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,3),4)
replicates_s4<-replicates_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,3,7,19),2)
replicates_s8<-replicates_analysis(Nondyn_simulate=Nondyn_simulate_0.02,var_empirical_bio=Nondyn_bio_V0.01_B5T5_rep1000_0.02,perturb_raw=c(1,3,7,19,9,21,25,13),1)
########################################################################
Iden_plot<-data.frame(Perturb=rep(c("Rep1","Rep2","Rep4","Rep8"),each=93),Iden=c(replicates_s1,replicates_s2,replicates_s4,replicates_s8))
Iden_plot$Perturb<-factor(Iden_plot$Perturb, levels = c("Rep1","Rep2","Rep4","Rep8"))

pp <- ggplot(Iden_plot, aes(x=Perturb,y=Iden)) + geom_boxplot(mapping=aes(fill=Perturb),outlier.shape = NA,color="grey20",alpha = 0.7) 
pp+labs(x = "",y = "")+
  geom_hline(yintercept=0,color = adjustcolor( "orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
  geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),color = "dodgerblue4",alpha = 0.4)+
  scale_fill_manual(values=c("dodgerblue4", "dodgerblue4", "dodgerblue4", "dodgerblue4"))+
  theme_bw()+
  theme(axis.text.x=element_text(size=rel(0),face="bold"),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Iden_replicates8_revise.pdf",width=5,height=4)
