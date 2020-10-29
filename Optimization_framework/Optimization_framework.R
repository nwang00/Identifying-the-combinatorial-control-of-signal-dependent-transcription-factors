# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Optimization_framework")
rm(list=ls())
library("MASS")
library(deSolve)
library(signal)

source("./ODE_functions.R")
Para_file_path<- commandArgs(trailingOnly=TRUE)
Paras_file<-read.table(file=Para_file_path)
#####################################################
##########################Process TF Activities###########################
#####################################################
TF_time_ca<-as.numeric(unlist(strsplit(Paras_file[4,2],",")))
TF_Activity_file<-read.table(file=Paras_file[1,2],header = T)
Num_pertb<-as.numeric(Paras_file[6,2])
int_length<-10^-3
Basal_level<-0

TF1_act<-vector()
TF2_act<-vector()
TF3_act<-vector()
for(TF_act_ind in 1:Num_pertb){
  TF1_act<-cbind(TF1_act,pchip(TF_time_ca, TF_Activity_file[,TF_act_ind], seq(0,(max(TF_time_ca)+20),int_length)))
  TF2_act<-cbind(TF2_act,pchip(TF_time_ca, TF_Activity_file[,(TF_act_ind+Num_pertb)], seq(0,(max(TF_time_ca)+20),int_length)))
  TF3_act<-cbind(TF3_act,pchip(TF_time_ca, TF_Activity_file[,(TF_act_ind+2*Num_pertb)], seq(0,(max(TF_time_ca)+20),int_length)))
}

Basal_level_TF<-0.05
TF1_act[TF1_act<Basal_level_TF]<-0
TF2_act[TF2_act<Basal_level_TF]<-0
TF3_act[TF3_act<Basal_level_TF]<-0
#####################################################
##########################Process Gene Exp###########################
#####################################################
Gene_exp_file<-read.table(file=Paras_file[2,2],header = T)
Gene_exp_time<-as.numeric(unlist(strsplit(Paras_file[5,2],",")))
exp_data<-as.numeric(Gene_exp_file)

##add sudo counts for zero counts
exp_data[exp_data==0]<-0.1

exp_data_WT_max<-max(exp_data[1:length(Gene_exp_time)])
num_sim_point<-2*Num_pertb*(length(Gene_exp_time)-1)+length(Gene_exp_time)
exp_data_basal<-rep(exp_data[seq(1,(Num_pertb*length(Gene_exp_time)),length(Gene_exp_time))],each=num_sim_point)
####################################################
#############Error for RNA
####################################################
# load("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Exp_fit/WT_LipidA_GeneExp.Prior.RData")
# load("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Exp_fit/MYD88_LipidA_GeneExp.Prior.RData")
# load("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Exp_fit/MAPKi_LipidA_GeneExp.Prior.RData")
# 
# var_est_g_list<-rbind(WT_LipidA_GeneExp.Prior$var_est_g,WT_LipidA_GeneExp.Prior$var_est_g,MYD88_LipidA_GeneExp.Prior$var_est_g,MAPKi_LipidA_GeneExp.Prior$var_est_g,WT_LipidA_GeneExp.Prior$var_est_g)
# var_est_n_list<-rbind(WT_LipidA_GeneExp.Prior$var_est_n,WT_LipidA_GeneExp.Prior$var_est_n,MYD88_LipidA_GeneExp.Prior$var_est_n,MAPKi_LipidA_GeneExp.Prior$var_est_n,WT_LipidA_GeneExp.Prior$var_est_n)

Gene_uncertainty_file<-t(read.table(file=Paras_file[3,2],header = T)[,-1])
var_est_g_list<-Gene_uncertainty_file[,1]
var_est_n_list<-Gene_uncertainty_file[,-1]
###############################calculate slope#################
slope_ind<-round(2*(var_est_n_list[,4])^0.5)
slope_ind[slope_ind>max(Gene_exp_time)]<-max(Gene_exp_time)
slope_ind[slope_ind<1]<-1

caRNA_time_org<-c(Gene_exp_time,as.vector(t(rbind(rep(Gene_exp_time[-1],each=2),rep(Gene_exp_time[-1],each=2))[rep(1,Num_pertb),]+c(-slope_ind,slope_ind))))
caRNA_time_org[caRNA_time_org<0]<-0
caRNA_time_org[caRNA_time_org>(max(Gene_exp_time)+20)]<-(max(Gene_exp_time)+20)
slope_interval_n<-(caRNA_time_org[seq((length(Gene_exp_time)+2),length(caRNA_time_org),2)]-caRNA_time_org[seq((length(Gene_exp_time)+1),length(caRNA_time_org),2)])

caRNA_time<-sort(caRNA_time_org)
back_order<-rep(order(order(caRNA_time_org)),Num_pertb)+rep(seq(0,num_sim_point*(Num_pertb-1),num_sim_point),each=num_sim_point)
################################################
################################################
data_ind<-rep(c(1:length(Gene_exp_time)),Num_pertb)+rep(seq(0,num_sim_point*(Num_pertb-1),num_sim_point),each=length(Gene_exp_time))
slope_list_n<-rep(0,length(Gene_exp_time)*Num_pertb)
slop_list_ind<-c(1:(length(Gene_exp_time)*Num_pertb))[-seq(1,length(Gene_exp_time)*Num_pertb,length(Gene_exp_time))]

slop_cal_ind1_n<-(c((length(Gene_exp_time)+1):num_sim_point)+rep((c(0:(Num_pertb-1))*num_sim_point),each=(2*(length(Gene_exp_time)-1))))[seq(1,(Num_pertb*(length(Gene_exp_time)-1)*2),2)]
slop_cal_ind2_n<-(c((length(Gene_exp_time)+1):num_sim_point)+rep((c(0:(Num_pertb-1))*num_sim_point),each=(2*(length(Gene_exp_time)-1))))[seq(2,(Num_pertb*(length(Gene_exp_time)-1)*2),2)]

var_est_g_cal<-var_est_g_list[rep(c(1:Num_pertb),each=length(Gene_exp_time))]
var_est_n_cal<-var_est_n_list[rep(c(1:Num_pertb),each=length(Gene_exp_time)),]

ind_fun<-function(x){
  return((x<=3))
}
##############################
# model_data<-rep(caRNA_time,5)

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
sample_size<-100

K1p_raw<-replicate(8,runif(sample_size_full, min=-2, max=5))
K2p_raw<-replicate(8,runif(sample_size_full, min=-2, max=5))
K3p_raw<-replicate(8,runif(sample_size_full, min=-2, max=5))
K5p<-replicate(8,runif(sample_size, min=-2, max=1))
#####################################################
#########Remove parameters combinations for weak gate output####
#####################################################
weak_threshold<-0.166

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
                                 times = Gene_exp_time,func =get(paste("P1_GRN",GRN_ind,"_fit",sep="")),parms = c(Kd1=10^(logKd1),Kd2=10^(logKd2),Kd3=10^(logKd3),ksyn=10^(logksyn),kdeg=10^(logkdeg)),hmax=1)[,2]))
}

for(GRN_ind in 1:8){
  tmp_opt<-vector()
  for(i in 1:sample_size) {
    Ksyn_sample<-ksyn_cal(logKd1 = K1p[i,GRN_ind],logKd2 = K2p[i,GRN_ind],logKd3 = K3p[i,GRN_ind],logksyn=0,logkdeg=K5p[i],Basal_level = Basal_level,exp_data_WT_max=exp_data_WT_max)
    tmp_opt<-c(tmp_opt,optim(c(K1p[i,GRN_ind],K2p[i,GRN_ind],K3p[i,GRN_ind],log10(Ksyn_sample),K5p[i]),method="L-BFGS-B",lower = c(rep(-2,3),-5,-2),upper = c(rep(5,3),6,1),P1_fun_GRNs,hessian=F))
  } 
  assign(paste("dyn_Gene_Logic",GRN_ind,sep=""),tmp_opt)
  save(list=paste("dyn_Gene_Logic",GRN_ind,sep=""),file = paste("dyn_Gene_Logic",GRN_ind,".RData",sep=""))
}
########################################################################
################################Examine results########################################
########################################################################
mse_GRNs<-vector()
par_GRNs<-vector()

for(GRN_ind in 1:8){
      sim_res_GRS<-get(paste("dyn_Gene_Logic",GRN_ind,sep=""))
      min_ind<-which(unlist(sim_res_GRS[seq(2,(5*sample_size),5)])==min(unlist(sim_res_GRS[seq(2,(5*sample_size),5)])))[1]
      mse_GRNs<-rbind(mse_GRNs,unlist(sim_res_GRS[min_ind*5-3]))
      par_GRNs[length(par_GRNs)+1]<-list(unlist(sim_res_GRS[min_ind*5-4]))
}
########################################################################
###############################Print results#########################################
########################################################################
result_print<-cbind(mse_GRNs,t(matrix(10^(unlist(par_GRNs)),nrow=5)))

rownames(result_print)<-c("TF1 AND TF2 AND TF3",
                        "TF1 OR TF2 OR TF3",
                        "TF1 OR TF2 AND TF3",
                        "TF1 AND (TF2 OR TF3)",
                        "TF2 OR TF1 AND TF3",
                        "TF2 AND (TF1 OR TF3)",
                        "TF3 OR TF1 AND TF2",
                        "TF3 AND (TF1 OR TF2)")
colnames(result_print)<-c("Negative_Log_Likelihood","Kd1(TF1)","Kd2(TF2)","Kd3(TF3)","ksyn","kdeg")

print("Fitting results of logic gates (ordered from good to bad fit)")
print(result_print[order(result_print[,1]),])
#######################Converted matrix#########################################
RS_matrix<-rbind(c(NA,NA,NA),c(NA,NA,NA))[rep(1,8),]
RS_raw<-result_print[,2:4]
RS_matrix[which(RS_raw<=0.5)]<-"Strong"
RS_matrix[intersect(which(RS_raw<=5),which(RS_raw>0.5))]<-"Medium"
RS_matrix[intersect(which(RS_raw<=100),which(RS_raw>5))]<-"Weak"
RS_matrix[which(RS_raw>100)]<-"Null"

result_convert<-as.data.frame(result_print)
result_convert[,2:4]<-RS_matrix
print("Results with converted regulation strength")
print(result_convert[order(result_convert[,1]),])


write.table(result_print[order(result_print[,1]),],file = "Optimization_output.txt",row.names = c("TF1 AND TF2 AND TF3",
                                                                                                  "TF1 OR TF2 OR TF3",
                                                                                                  "TF1 OR TF2 AND TF3",
                                                                                                  "TF1 AND (TF2 OR TF3)",
                                                                                                  "TF2 OR TF1 AND TF3",
                                                                                                  "TF2 AND (TF1 OR TF3)",
                                                                                                  "TF3 OR TF1 AND TF2",
                                                                                                  "TF3 AND (TF1 OR TF2)"),
            quote = F,col.names =c("Negative_Log_Likelihood","Kd1(TF1)","Kd2(TF2)","Kd3(TF3)","ksyn","kdeg"))
