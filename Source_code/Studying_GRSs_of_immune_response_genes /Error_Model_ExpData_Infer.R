# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Exp_test")
rm(list=ls())
# library("mixtools")
library(signal)
library(MASS)
# # library(ggplot2)
# library(viridis)
# library(grid)
########################################################################
####################################function
########################################################################
NLL_Prior_g<-function(data_point,mean_est,X,gamma_ind_all){
  NLL_gamma<-(-sum(log(dgamma(data_point[gamma_ind_all],shape=(mean_est[gamma_ind_all])^2/(10^(X[1])),scale = (10^(X[1]))/mean_est[gamma_ind_all])+10^(-20))))
  if(is.infinite((NLL_gamma))){
    return(10^8)
  }else{
    return((NLL_gamma))
  }
}
########################################################################
NLL_Prior_TV_n<-function(data_point,mean_est,slop_square,X,normal_ind_all){
  slope_ind<-round(2*((2*10^(X[4]))^0.5))
  if(slope_ind>60) slope_ind<-60
  if(slope_ind<1) slope_ind<-1
  
  f_mean<-10^(X[1])+10^(X[2])*mean_est+10^(X[3])*(mean_est)^2
  var_total<-f_mean+slop_square[slope_ind,]*10^(X[4])
  
  NLL_normal<-(-sum(log(dnorm(data_point[normal_ind_all],mean=mean_est[normal_ind_all],sd=sqrt(var_total[normal_ind_all]))+10^(-20))))
  
  if(is.infinite((NLL_normal))){
    return(10^8)
  }else{
    return((NLL_normal))
  }
}
########################################################################
########################Prior estimation#########################
############################Teporal Value###############################
########################################################################
Prior_TV_estimate<-function(data.raw,caRNA_Slope_Full_list){
  # data.raw<-WT_LipidA_GeneExp
  # caRNA_Slope_Full_list<-caRNA_Slope_bio_V0.01_B5T5_rep1000
  data.mean<-apply(data.raw,1,mean)
  
  normal_ind<-which(data.mean>3)
  zero_ind<-which(apply(data.raw,1,min)==0)
  gamma_ind<-setdiff(setdiff(c(1:length(data.mean)),normal_ind),zero_ind)
  
  gamma_ind_all<-c(gamma_ind,(gamma_ind+dim(data.raw)[1]))
  normal_ind_all<-c(normal_ind,(normal_ind+dim(data.raw)[1]))
  
  K1p<-runif(1000, min=-20, max=10)
  K2p<-runif(1000, min=-20, max=10)
  K3p<-runif(1000, min=-20, max=10)
  K4p<-runif(1000, min=-20, max=10)
  
  RNA_gw_list<-vector()
  print("Start estimation normal")
  pb = txtProgressBar(min = 1, max = 200)
  
  for(i in 1:200){
    RNA_LH<-optim(c(K1p[i],K2p[i],K3p[i],K4p[i]),method="L-BFGS-B",fn=NLL_Prior_TV_n,lower=rep((-20),4),upper=rep((10),4),normal_ind_all=normal_ind_all,
                  data_point=c(data.raw[,1],data.raw[,2]),slop_square=cbind((caRNA_Slope_Full_list)^2,(caRNA_Slope_Full_list)^2),
                  mean_est=rep(data.mean,2))
    
    RNA_gw_list<-c(RNA_gw_list,RNA_LH)
    setTxtProgressBar(pb,i)
  }
  
  var_est_n<-2*10^(RNA_gw_list[order(unlist(RNA_gw_list[seq(2,5*1000,5)]))[1]*5-4]$par)
  ######slope estimate#########
  slope_ind<-round(2*(var_est_n[4])^0.5)
  if(slope_ind>60) slope_ind<-60
  if(slope_ind<1) slope_ind<-1
  
  f_mean<-var_est_n[1]+var_est_n[2]*data.mean+var_est_n[3]*(data.mean)^2
  var_est_each_n<-f_mean+(caRNA_Slope_Full_list[slope_ind,])^2*var_est_n[4]
  ###############################################################
  ###############################################################
  ###############################################################
  K1p<-runif(1000, min=-20, max=10)
  
  RNA_gw_list<-vector()
  print("Start estimation gamma")
  pb = txtProgressBar(min = 1, max = 200)
  
  for(i in 1:200){
    RNA_LH<-optim(c(K1p[i]),method="L-BFGS-B",fn=NLL_Prior_g,lower=(-20),upper=(10),gamma_ind_all=gamma_ind_all,
                  data_point=c(data.raw[,1],data.raw[,2]),
                  mean_est=rep(data.mean,2))
    
    RNA_gw_list<-c(RNA_gw_list,RNA_LH)
    setTxtProgressBar(pb,i)
  }
  
  var_est_g<-2*10^(RNA_gw_list[order(unlist(RNA_gw_list[seq(2,5*1000,5)]))[1]*5-4]$par)
  ###############################################################
  var_est_each<-(c(1:length(data.mean))%in%gamma_ind)*var_est_g+(c(1:length(data.mean))%in%normal_ind)*var_est_each_n
  ###############################################################
  return(list(data.mean=data.mean,gamma_ind=gamma_ind,normal_ind=normal_ind,var_est_g=var_est_g,var_est_n=var_est_n,var_est_each=var_est_each))
}
#################################Slope Estimate######################
Slope_estimate_f<-function(caRNA_raw.mean){
  caRNA_Slope_Full_list<-vector()
  for(slope_delta in 1:60){
    caRNA_Slope_list<-rep(0,length(caRNA_raw.mean))
    for(i in 1:(length(caRNA_raw.mean)/4)){
      ind_delta<-c((15-slope_delta),(15+slope_delta),(30-slope_delta),(30+slope_delta),(60-slope_delta),(60+slope_delta))
      ind_delta[ind_delta<0]<-0
      ind_delta[ind_delta>80]<-80
      temp<-pchip(c(0,15,30,60),caRNA_raw.mean[(i*4-3):(i*4)],ind_delta)
      
      caRNA_Slope_list[(i*4-2):(i*4)]<-c((temp[2]-temp[1]),(temp[4]-temp[3]),(temp[6]-temp[5]))/(c((ind_delta[2]-ind_delta[1]),(ind_delta[4]-ind_delta[3]),(ind_delta[6]-ind_delta[5])))
      # slope for non activated to be 0
      if(all(caRNA_raw.mean[(i*4-3):(i*4)]<3)){
        caRNA_Slope_list[(i*4-2):(i*4)]<-0
      }
    }
    # #set negative slope to be 0
    caRNA_Slope_list[caRNA_Slope_list<0]<-0
    caRNA_Slope_Full_list<-rbind(caRNA_Slope_Full_list,caRNA_Slope_list)
  }
  return(caRNA_Slope_Full_list)
}
########################################################################
#############################Load Data V4#######################################
################################################################################
load("RPKM_Smale.RData")
WT_LipidA_GeneExp<-cbind(matrix(t(RPKM_Smale[,29:32]),ncol=1),matrix(t(RPKM_Smale[,34:37]),ncol=1))
WT_LipidA_GeneExp.mean<-apply(WT_LipidA_GeneExp,1,mean)

MYD88_LipidA_GeneExp<-cbind(matrix(t(RPKM_Smale[,c(100,104,106,108)]),ncol=1),matrix(t(RPKM_Smale[,c(101,105,107,109)]),ncol=1))
MYD88_LipidA_GeneExp.mean<-apply(MYD88_LipidA_GeneExp,1,mean)
####################################################################
WT_lipidA_mean_matirx<-(RPKM_Smale[,29:32]+RPKM_Smale[,34:37])/2
Control_MAPKi_lipidA_mean_matirx<-(RPKM_Smale[,c(71,75,77,81)]+RPKM_Smale[,c(70,74,76,80)])/2
norm_factor_raw<-WT_lipidA_mean_matirx/(Control_MAPKi_lipidA_mean_matirx+0.1)
norm_factor<-norm_factor_raw
norm_factor[norm_factor==0]<-1
norm_factor_final<-norm_factor[,c(1,1,2,2,3,3,4,4)]

MAPKi_LipidA_GeneExp_raw<-RPKM_Smale[,c(40,39,44,43,46,45,50,49)]*norm_factor_final
MAPKi_LipidA_GeneExp<-cbind(matrix(t(MAPKi_LipidA_GeneExp_raw[,c(1,3,5,7)]),ncol=1),matrix(t(MAPKi_LipidA_GeneExp_raw[,c(2,4,6,8)]),ncol=1))
MAPKi_LipidA_GeneExp.mean<-apply(MAPKi_LipidA_GeneExp,1,mean)
####################################################################
###########Slope Estimate and Plot
####################################################################
WT_LipidA_GeneExp.slope<-Slope_estimate_f(WT_LipidA_GeneExp.mean)
MYD88_LipidA_GeneExp.slope<-Slope_estimate_f(MYD88_LipidA_GeneExp.mean)
MAPKi_LipidA_GeneExp.slope<-Slope_estimate_f(MAPKi_LipidA_GeneExp.mean)
save.image("Slope_exp.RData")
####################################################################
#######################################Variance Estimation############################
#####################################all parameters###############################
WT_LipidA_GeneExp.Prior<-Prior_TV_estimate(WT_LipidA_GeneExp,WT_LipidA_GeneExp.slope)
save(WT_LipidA_GeneExp.Prior,file = "WT_LipidA_GeneExp.Prior.RData")

MYD88_LipidA_GeneExp.Prior<-Prior_TV_estimate(MYD88_LipidA_GeneExp,MYD88_LipidA_GeneExp.slope)
save(MYD88_LipidA_GeneExp.Prior,file = "MYD88_LipidA_GeneExp.Prior.RData")

MAPKi_LipidA_GeneExp.Prior<-Prior_TV_estimate(MAPKi_LipidA_GeneExp,MAPKi_LipidA_GeneExp.slope)
save(MAPKi_LipidA_GeneExp.Prior,file = "MAPKi_LipidA_GeneExp.Prior.RData")