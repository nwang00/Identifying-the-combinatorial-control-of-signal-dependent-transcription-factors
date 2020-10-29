# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics")
rm(list=ls())
# library("mixtools")
library(signal)
library(MASS)
library(ggplot2)
library(viridis)
library(grid)
########################################################################
########################################################################
####################################function
########################################################################
########################################################################
#######################################Basal Estimate#################################
NLL_MLE_g<-function(data_point,mean_est,X){
  var_total<-10^(X[1])
  NLL<-(-sum(log(dgamma(data_point, shape=(mean_est)^2/(var_total),scale = var_total/mean_est)+10^(-20))))
  if(is.infinite(NLL)){
    return(10^8)
  }else{
    return(NLL)
  }
}
######################################MLE:Value Temporal##################################
NLL_MLE_TV_n<-function(data_point,mean_est,slop_square,X){
  slope_ind<-round(2*((2*10^(X[4]))^0.5))
  if(slope_ind>60) slope_ind<-60
  if(slope_ind<1) slope_ind<-1
  f_mean<-10^(X[1])+10^(X[2])*mean_est+10^(X[3])*(mean_est)^2
  var_total<-f_mean+slop_square[slope_ind]*10^(X[4])
  
  NLL<-(-sum(log(dnorm(data_point,mean=mean_est,sd=sqrt(var_total))+10^(-20))))
  if(is.infinite(NLL)){
    return(10^8)
  }else{
    return(NLL)
  }
}
######################################MLE:Conventional##################################
NLL_MLE_Con_n<-function(data_point,mean_est,X){
  f_mean<-10^(X[1])+10^(X[2])*mean_est+10^(X[3])*(mean_est)^2
  
  NLL<-(-sum(log(dnorm(data_point,mean=mean_est,sd=sqrt(f_mean))+10^(-20))))
  if(is.infinite(NLL)){
    return(10^8)
  }else{
    return(NLL)
  }
}
########################################################################
########################################################################
########################################################################
NLL_Prior_g<-function(data_point,mean_est,X,gamma_ind_all){
  NLL_gamma<-(-sum(log(dgamma(data_point[gamma_ind_all],shape=(mean_est[gamma_ind_all])^2/(10^(X[1])),scale = (10^(X[1]))/mean_est[gamma_ind_all])+10^(-20))))
  # print((NLL_gamma+NLL_normal))
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
  
  # print((NLL_gamma+NLL_normal))
  if(is.infinite((NLL_normal))){
    return(10^8)
  }else{
    return((NLL_normal))
  }
}
########################################################################
NLL_Prior_Con_n<-function(data_point,mean_est,X,gamma_ind_all,normal_ind_all){
  f_mean<-10^(X[1])+10^(X[2])*mean_est+10^(X[3])*(mean_est)^2
  
  NLL_normal<-(-sum(log(dnorm(data_point[normal_ind_all],mean=mean_est[normal_ind_all],sd=sqrt(f_mean[normal_ind_all]))+10^(-20))))
  # print((NLL_gamma+NLL_normal))
  if(is.infinite((NLL_normal))){
    return(10^8)
  }else{
    return((NLL_normal))
  }
}
##################################MAP basal##################################
NLL_genewise_MAP_g<-function(data_point,mean_est,var_trend_g,sd_g,X){
  
  var_total<-10^(X[1])
  
  NLL<-(-sum(log(dgamma(data_point, shape=(mean_est)^2/(var_total),scale = (var_total)/mean_est)+10^(-20))))-sum(log(dnorm(X[1],mean=log10(var_trend_g[1]/2),sd=sd_g[1])+10^(-20)))
  if(is.infinite(NLL)){
    return(10^8)
  }else{
    return(NLL)
  }
}
##################################
NLL_genewise_MAP_TV_n<-function(data_point,mean_est,slop_square,var_trend_n,sd_n,X){
  # X<-RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-4]$par
  # X<-log10(var_trend/2)
  slope_ind<-round(2*((2*10^(X[4]))^0.5))
  if(slope_ind>60) slope_ind<-60
  if(slope_ind<1) slope_ind<-1
  
  f_mean<-10^(X[1])+10^(X[2])*mean_est+10^(X[3])*(mean_est)^2
  var_total<-f_mean+slop_square[slope_ind]*10^(X[4])
  
  NLL<-(-sum(log(dnorm(data_point,mean=mean_est,sd=sqrt(var_total))+10^(-20))))-sum(log(dnorm(X[1],mean=log10(var_trend_n[1]/2),sd=sd_n[1])+10^(-20)))-sum(log(dnorm(X[2],mean=log10(var_trend_n[2]/2),sd=sd_n[2])+10^(-20)))-sum(log(dnorm(X[3],mean=log10(var_trend_n[3]/2),sd=sd_n[3])+10^(-20)))-sum(log(dnorm(X[4],mean=log10(var_trend_n[4]/2),sd=sd_n[4])+10^(-20)))
  # print(NLL)
  if(is.infinite(NLL)){
    return(10^8)
  }else{
    return(NLL)
  }
}
####################################################################
NLL_genewise_MAP_Con_n<-function(data_point,mean_est,var_trend_n,sd_n,X){
  
  f_mean<-10^(X[1])+10^(X[2])*mean_est+10^(X[3])*(mean_est)^2
  var_total<-f_mean
  
  NLL<-(-sum(log(dnorm(data_point,mean=mean_est,sd=sqrt(var_total))+10^(-20))))-sum(log(dnorm(X[1],mean=log10(var_trend_n[1]/2),sd=sd_n[1])+10^(-20)))-sum(log(dnorm(X[2],mean=log10(var_trend_n[2]/2),sd=sd_n[2])+10^(-20)))-sum(log(dnorm(X[3],mean=log10(var_trend_n[3]/2),sd=sd_n[3])+10^(-20)))
  # print(NLL)
  if(is.infinite(NLL)){
    return(10^8)
  }else{
    return(NLL)
  }
}
####################################################################
##density plot function
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
####################################################################
####################################MLE Function: Temporal Value################################
####################################################################
MLE_TV_estimate<-function(data.raw,caRNA_Slope_Full_list){
  # data.raw<-caRNA_raw_BT50_parV_rep1000
  # caRNA_Slope_Full_list<-caRNA_Slope_BT50_parV_rep1000
  ind_non_zero<-intersect(which(data.raw[,1]!=0),which(data.raw[,2]!=0))
  data.mean<-apply(data.raw,1,mean)
  
  normal_ind<-intersect(which(data.mean>3),ind_non_zero)
  gamma_ind<-setdiff(c(1:length(data.mean)),normal_ind)
  
  var_est<-cbind(rep(0,length(data.mean)),rep(0,length(data.mean)),rep(0,length(data.mean)),rep(0,length(data.mean)))
  var_est_each<-rep(0,length(data.mean))
  var_est_quality<-rep(0,length(data.mean))
  
  if(length(gamma_ind)!=0){
    print("Start gamma estimate")
    Gene_var_g<-vector()
    if(length(gamma_ind)>1){
      pb = txtProgressBar(min = 1, max = length(gamma_ind)) 
    }
    
    for(j in 1:length(gamma_ind)){
      K1p<-runif(100, min=-20, max=10)
      
      
      RNA_gw_g_list<-vector()
      for(i in 1:30){
        RNA_LH_g<-optim(c(K1p[i]),method="L-BFGS-B",fn=NLL_MLE_g,lower=c(-20),upper=c(10),
                        data_point=c(data.raw[gamma_ind[j],1],data.raw[gamma_ind[j],2]),
                        mean_est=rep(data.mean[gamma_ind[j]],2))
        RNA_gw_g_list<-c(RNA_gw_g_list,RNA_LH_g)
      }
      
      var_est[gamma_ind[j],]<-c(2*10^(RNA_gw_g_list[order(unlist(RNA_gw_g_list[seq(2,5*100,5)]))[1]*5-4]$par),0,0,0)
      var_est_quality[gamma_ind[j]]<-RNA_gw_g_list[order(unlist(RNA_gw_g_list[seq(2,5*100,5)]))[1]*5-3]$value
      
      var_est_each[gamma_ind[j]]<-var_est[gamma_ind[j],1]
      
      setTxtProgressBar(pb,j)
      # if(j%%length(gamma_ind)==0){
      #   print(j)
      # }
    }
  }
  
  print("Start normal estimate")
  Gene_var_n<-vector()
  pb = txtProgressBar(min = 1, max = length(normal_ind)) 
  for(j in 1:length(normal_ind)){
    K1p<-runif(100, min=-20, max=10)
    K2p<-runif(100, min=-20, max=10)
    K3p<-runif(100, min=-20, max=10)
    K4p<-runif(100, min=-20, max=10)
    
    RNA_gw_n_list<-vector()
    for(i in 1:30){
      RNA_LH_n<-optim(c(K1p[i],K2p[i],K3p[i],K4p[i]),method="L-BFGS-B",fn=NLL_MLE_TV_n,lower=rep((-20),4),upper=rep((10),4),
                      data_point=c(data.raw[normal_ind[j],1],data.raw[normal_ind[j],2]),slop_square=rep((caRNA_Slope_Full_list[,normal_ind[j]])^2,2),
                      mean_est=rep(data.mean[normal_ind[j]],2))
      RNA_gw_n_list<-c(RNA_gw_n_list,RNA_LH_n)
    }
    var_est[normal_ind[j],]<-2*10^(RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-4]$par)
    var_est_quality[normal_ind[j]]<-RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-3]$value
    
    ######slope estimate#########
    slope_ind<-round(2*(var_est[normal_ind[j],4])^0.5)
    if(slope_ind>60) slope_ind<-60
    if(slope_ind<1) slope_ind<-1
    
    f_mean<-var_est[normal_ind[j],1]+var_est[normal_ind[j],2]*data.mean[normal_ind[j]]+var_est[normal_ind[j],3]*(data.mean[normal_ind[j]])^2+var_est[normal_ind[j],4]*(data.mean[normal_ind[j]])^3
    var_est_each[normal_ind[j]]<-f_mean+(caRNA_Slope_Full_list[slope_ind,normal_ind[j]])^2*var_est[normal_ind[j],4]
    
    setTxtProgressBar(pb,j)
  }
  
  return(list(data.mean=data.mean,gamma_ind=gamma_ind,normal_ind=normal_ind,var_est=var_est,var_est_each=var_est_each,var_est_quality=var_est_quality))
}
#################Conventional
MLE_Con_estimate<-function(data.raw){
  # data.raw<-caRNA_raw_BT50_parV_rep1000
  ind_non_zero<-intersect(which(data.raw[,1]!=0),which(data.raw[,2]!=0))
  data.mean<-apply(data.raw,1,mean)
  
  normal_ind<-intersect(which(data.mean>3),ind_non_zero)
  gamma_ind<-setdiff(c(1:length(data.mean)),normal_ind)
  
  var_est<-cbind(rep(0,length(data.mean)),rep(0,length(data.mean)),rep(0,length(data.mean)))
  var_est_each<-rep(0,length(data.mean))
  var_est_quality<-rep(0,length(data.mean))
  
  if(length(gamma_ind)!=0){
    print("Start gamma estimate")
    Gene_var_g<-vector()
    if(length(gamma_ind)>1){
      pb = txtProgressBar(min = 1, max = length(gamma_ind)) 
    }
    
    for(j in 1:length(gamma_ind)){
      K1p<-runif(100, min=-20, max=10)
      
      RNA_gw_g_list<-vector()
      for(i in 1:30){
        RNA_LH_g<-optim(c(K1p[i]),method="L-BFGS-B",fn=NLL_MLE_g,lower=c(-20),upper=c(10),
                        data_point=c(data.raw[gamma_ind[j],1],data.raw[gamma_ind[j],2]),
                        mean_est=rep(data.mean[gamma_ind[j]],2))
        RNA_gw_g_list<-c(RNA_gw_g_list,RNA_LH_g)
      }
      
      var_est[gamma_ind[j],]<-c(2*10^(RNA_gw_g_list[order(unlist(RNA_gw_g_list[seq(2,5*100,5)]))[1]*5-4]$par),0,0)
      var_est_quality[gamma_ind[j]]<-RNA_gw_g_list[order(unlist(RNA_gw_g_list[seq(2,5*100,5)]))[1]*5-3]$value
      
      var_est_each[gamma_ind[j]]<-var_est[gamma_ind[j],1]
      
      setTxtProgressBar(pb,j)
      # if(j%%length(gamma_ind)==0){
      #   print(j)
      # }
      # print(j)
    }
  }
  
  print("Start normal estimate")
  Gene_var_n<-vector()
  pb = txtProgressBar(min = 1, max = length(normal_ind)) 
  for(j in 1:length(normal_ind)){
    K1p<-runif(100, min=-20, max=10)
    K2p<-runif(100, min=-20, max=10)
    K3p<-runif(100, min=-20, max=10)
    
    RNA_gw_n_list<-vector()
    for(i in 1:30){
      RNA_LH_n<-optim(c(K1p[i],K2p[i],K3p[i]),method="L-BFGS-B",fn=NLL_MLE_Con_n,lower=rep((-20),3),upper=rep((10),3),
                      data_point=c(data.raw[normal_ind[j],1],data.raw[normal_ind[j],2]),
                      mean_est=rep(data.mean[normal_ind[j]],2))
      RNA_gw_n_list<-c(RNA_gw_n_list,RNA_LH_n)
    }
    
    var_est[normal_ind[j],]<-2*10^(RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-4]$par)
    var_est_quality[normal_ind[j]]<-RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-3]$value
    
    f_mean<-var_est[normal_ind[j],1]+var_est[normal_ind[j],2]*data.mean[normal_ind[j]]+var_est[normal_ind[j],3]*(data.mean[normal_ind[j]])^2
    var_est_each[normal_ind[j]]<-f_mean
    
    setTxtProgressBar(pb,j)
  }
  return(list(data.mean=data.mean,gamma_ind=gamma_ind,normal_ind=normal_ind,var_est=var_est,var_est_each=var_est_each,var_est_quality=var_est_quality))
}
########################################################################
########################################################################
########################Prior estimation#########################
############################Teporal Value###############################
########################################################################
Prior_TV_estimate<-function(data.raw,caRNA_Slope_Full_list){
  # data.raw<-caRNA_raw_bio_V0.01_BT50_rep1000
  # caRNA_Slope_Full_list<-caRNA_Slope_bio_V0.01_BT50_rep1000
  data.mean<-apply(data.raw,1,mean)
  
  normal_ind<-which(data.mean>3)
  gamma_ind<-setdiff(c(1:length(data.mean)),normal_ind)
  
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
############################Conventional ###############################
Prior_Con_estimate<-function(data.raw){
  # data.raw<-caRNA_raw_BT50_parV_rep1000
  data.mean<-apply(data.raw,1,mean)
  
  normal_ind<-which(data.mean>3)
  gamma_ind<-setdiff(c(1:length(data.mean)),normal_ind) 
  gamma_ind_all<-c(gamma_ind,(gamma_ind+dim(data.raw)[1]))
  normal_ind_all<-c(normal_ind,(normal_ind+dim(data.raw)[1]))
  
  K1p<-runif(1000, min=-20, max=10)
  K2p<-runif(1000, min=-20, max=10)
  K3p<-runif(1000, min=-20, max=10)
  
  RNA_gw_list<-vector()
  print("Start estimation")
  pb = txtProgressBar(min = 1, max = 200)
  for(i in 1:200){
    RNA_LH<-optim(c(K1p[i],K2p[i],K3p[i]),method="L-BFGS-B",fn=NLL_Prior_Con_n,lower=rep((-20),3),upper=rep((10),3),normal_ind_all=normal_ind_all,
                  data_point=c(data.raw[,1],data.raw[,2]),
                  mean_est=rep(data.mean,2))
    
    RNA_gw_list<-c(RNA_gw_list,RNA_LH)
    setTxtProgressBar(pb,i)
  }
  
  var_est_n<-2*10^(RNA_gw_list[order(unlist(RNA_gw_list[seq(2,5*1000,5)]))[1]*5-4]$par)
  
  f_mean<-var_est_n[1]+var_est_n[2]*data.mean+var_est_n[3]*(data.mean)^2
  var_est_each_n<-f_mean
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
########################################################################
########################################################################
########################MAP estimation#########################
############################Teporal Value###############################
MAP_TV_estimate<-function(data.raw,caRNA_Slope_Full_list,var_trend_g,var_trend_n,sd_coeff){
  # data.raw<-caRNA_raw_bio_V0.01_BT50_rep1000
  # caRNA_Slope_Full_list<-caRNA_Slope_bio_V0.01_BT50_rep1000
  # var_trend_n=Prior_TV_bio_V0.01_BT50_rep1000_M0.75_imp$var_est_n
  # var_trend_g=Prior_TV_bio_V0.01_BT50_rep1000_M0.75_imp$var_est_g
  # sd_coeff<-M0.75
  data.mean<-apply(data.raw,1,mean)
  
  normal_ind<-which(data.mean>3)
  gamma_ind<-setdiff(c(1:length(data.mean)),normal_ind)
  
  gamma_ind_all<-c(gamma_ind,(gamma_ind+dim(data.raw)[1]))
  normal_ind_all<-c(normal_ind,(normal_ind+dim(data.raw)[1]))
  
  sd_prior_g<-c(sd_coeff*abs(log10(var_trend_g[1]/2)))
  sd_prior_n<-c(sd_coeff*abs(log10(var_trend_n[1]/2)),sd_coeff*abs(log10(var_trend_n[2]/2)),sd_coeff*abs(log10(var_trend_n[3]/2)),sd_coeff*abs(log10(var_trend_n[4]/2)))
  
  var_est<-cbind(rep(0,length(data.mean)),rep(0,length(data.mean)),rep(0,length(data.mean)),rep(0,length(data.mean)))
  var_est_each<-rep(0,length(data.mean))
  var_est_quality<-rep(0,length(data.mean))
  
  if(length(gamma_ind)!=0){
    print("Start gamma estimate")
    Gene_var_g<-vector()
    if(length(gamma_ind)>1){
      pb = txtProgressBar(min = 1, max = length(gamma_ind)) 
    }
    
    for(j in 1:length(gamma_ind)){
      K1p<-rnorm(100,mean=log10(var_trend_g[1]/2),sd=sd_prior_g[1])
      
      RNA_gw_g_list<-vector()
      for(i in 1:30){
        RNA_LH_g<-optim(c(K1p[i]),method="L-BFGS-B",fn=NLL_genewise_MAP_g,lower=c(-20),upper=c(10),
                        data_point=c(data.raw[gamma_ind[j],1],data.raw[gamma_ind[j],2]),
                        mean_est=rep(data.mean[gamma_ind[j]],2),var_trend_g=var_trend_g,sd_g=sd_prior_g)
        
        
        RNA_gw_g_list<-c(RNA_gw_g_list,RNA_LH_g)
      }
      
      var_est[gamma_ind[j],]<-c(2*10^(RNA_gw_g_list[order(unlist(RNA_gw_g_list[seq(2,5*100,5)]))[1]*5-4]$par),0,0,0)
      var_est_quality[gamma_ind[j]]<-RNA_gw_g_list[order(unlist(RNA_gw_g_list[seq(2,5*100,5)]))[1]*5-3]$value
      ######slope estimate#########
      var_est_each[gamma_ind[j]]<-var_est[gamma_ind[j],]
      setTxtProgressBar(pb,j)
      # if(j%%length(gamma_ind)==0){
      #   print(j)
      # }
    }
  }
  
  print("Start normal estimate")
  Gene_var_n<-vector()
  pb = txtProgressBar(min = 1, max = length(normal_ind)) 
  for(j in 1:length(normal_ind)){
    K1p<-rnorm(100,mean=log10(var_trend_n[1]/2),sd=sd_prior_n[1])
    K2p<-rnorm(100,mean=log10(var_trend_n[2]/2),sd=sd_prior_n[2])
    K3p<-rnorm(100,mean=log10(var_trend_n[3]/2),sd=sd_prior_n[3])
    K4p<-rnorm(100,mean=log10(var_trend_n[4]/2),sd=sd_prior_n[4])
    
    RNA_gw_n_list<-vector()
    for(i in 1:30){
      RNA_LH_n<-optim(c(K1p[i],K2p[i],K3p[i],K4p[i]),method="L-BFGS-B",fn=NLL_genewise_MAP_TV_n,lower=rep((-20),4),upper=rep((10),4),
                      data_point=c(data.raw[normal_ind[j],1],data.raw[normal_ind[j],2]),slop_square=rep((caRNA_Slope_Full_list[,normal_ind[j]])^2,2),
                      mean_est=rep(data.mean[normal_ind[j]],2),var_trend_n=var_trend_n,sd_n=sd_prior_n)
      RNA_gw_n_list<-c(RNA_gw_n_list,RNA_LH_n)
    }
    var_est[normal_ind[j],]<-2*10^(RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-4]$par)
    var_est_quality[normal_ind[j]]<-RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-3]$value
    ######slope estimate#########
    slope_ind<-round(2*(var_est[normal_ind[j],4])^0.5)
    if(slope_ind>60) slope_ind<-60
    if(slope_ind<1) slope_ind<-1
    
    f_mean<-var_est[normal_ind[j],1]+var_est[normal_ind[j],2]*data.mean[normal_ind[j]]+var_est[normal_ind[j],3]*(data.mean[normal_ind[j]])^2
    var_est_each[normal_ind[j]]<-f_mean+(caRNA_Slope_Full_list[slope_ind,normal_ind[j]])^2*var_est[normal_ind[j],4]
    
    setTxtProgressBar(pb,j)
  }
  
  return(list(data.mean=data.mean,gamma_ind=gamma_ind,normal_ind=normal_ind,var_est=var_est,var_est_each=var_est_each,var_est_quality=var_est_quality))
}
############################Conventional###############################
MAP_Con_estimate<-function(data.raw,var_trend_g,var_trend_n,sd_coeff){
  # data.raw<-caRNA_raw_BT50_parV_rep1000
  # 
  # var_trend=Prior_Con_BT50_parV_rep1000$var_est
  # sd_coeff<-0.0001
  data.mean<-apply(data.raw,1,mean)
  normal_ind<-which(data.mean>3)
  gamma_ind<-setdiff(c(1:length(data.mean)),normal_ind)
  
  gamma_ind_all<-c(gamma_ind,(gamma_ind+dim(data.raw)[1]))
  normal_ind_all<-c(normal_ind,(normal_ind+dim(data.raw)[1]))
  sd_prior_g<-c(sd_coeff*abs(log10(var_trend_g[1]/2)))
  sd_prior_n<-c(sd_coeff*abs(log10(var_trend_n[1]/2)),sd_coeff*abs(log10(var_trend_n[2]/2)),sd_coeff*abs(log10(var_trend_n[3]/2)))
  
  var_est<-cbind(rep(0,length(data.mean)),rep(0,length(data.mean)),rep(0,length(data.mean)))
  var_est_each<-rep(0,length(data.mean))
  var_est_quality<-rep(0,length(data.mean))
  
  if(length(gamma_ind)!=0){
    print("Start gamma estimate")
    Gene_var_g<-vector()
    if(length(gamma_ind)>1){
      pb = txtProgressBar(min = 1, max = length(gamma_ind)) 
    }
    
    for(j in 1:length(gamma_ind)){
      K1p<-rnorm(100,mean=log10(var_trend_g[1]/2),sd=sd_prior_g[1])
      
      RNA_gw_g_list<-vector()
      for(i in 1:30){
        RNA_LH_g<-optim(c(K1p[i]),method="L-BFGS-B",fn=NLL_genewise_MAP_g,lower=(-20),upper=(10),
                        data_point=c(data.raw[gamma_ind[j],1],data.raw[gamma_ind[j],2]),
                        mean_est=rep(data.mean[gamma_ind[j]],2),var_trend_g=var_trend_g,sd_g=sd_prior_g)
        RNA_gw_g_list<-c(RNA_gw_g_list,RNA_LH_g)
      }
      
      var_est[gamma_ind[j],]<-c(2*10^(RNA_gw_g_list[order(unlist(RNA_gw_g_list[seq(2,5*100,5)]))[1]*5-4]$par),0,0)
      var_est_quality[gamma_ind[j]]<-RNA_gw_g_list[order(unlist(RNA_gw_g_list[seq(2,5*100,5)]))[1]*5-3]$value
      ######slope estimate#########
      var_est_each[gamma_ind[j]]<-var_est[gamma_ind[j],1]
      
      setTxtProgressBar(pb,j)
      # if(j%%length(gamma_ind)==0){
      #   print(j)
      # }
    }
  }
  
  print("Start normal estimate")
  Gene_var_n<-vector()
  pb = txtProgressBar(min = 1, max = length(normal_ind)) 
  for(j in 1:length(normal_ind)){
    K1p<-rnorm(100,mean=log10(var_trend_n[1]/2),sd=sd_prior_n[1])
    K2p<-rnorm(100,mean=log10(var_trend_n[2]/2),sd=sd_prior_n[2])
    K3p<-rnorm(100,mean=log10(var_trend_n[3]/2),sd=sd_prior_n[3])
    
    RNA_gw_n_list<-vector()
    for(i in 1:30){
      RNA_LH_n<-optim(c(K1p[i],K2p[i],K3p[i]),method="L-BFGS-B",fn=NLL_genewise_MAP_Con_n,lower=rep((-20),3),upper=rep((10),3),
                      data_point=c(data.raw[normal_ind[j],1],data.raw[normal_ind[j],2]),
                      mean_est=rep(data.mean[normal_ind[j]],2),var_trend_n=var_trend_n,sd_n=sd_prior_n)
      RNA_gw_n_list<-c(RNA_gw_n_list,RNA_LH_n)
    }
    var_est[normal_ind[j],]<-2*10^(RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-4]$par)
    var_est_quality[normal_ind[j]]<-RNA_gw_n_list[order(unlist(RNA_gw_n_list[seq(2,5*100,5)]))[1]*5-3]$value
    
    f_mean<-var_est[normal_ind[j],1]+var_est[normal_ind[j],2]*data.mean[normal_ind[j]]+var_est[normal_ind[j],3]*(data.mean[normal_ind[j]])^2
    var_est_each[normal_ind[j]]<-f_mean
    
    setTxtProgressBar(pb,j)
  }
  
  
  return(list(data.mean=data.mean,gamma_ind=gamma_ind,normal_ind=normal_ind,var_est=var_est,var_est_each=var_est_each,var_est_quality=var_est_quality))
}
#################################Slope Estimate######################
Slope_estimate_f<-function(caRNA_raw.mean){
  caRNA_Slope_Full_list<-vector()
  for(slope_delta in 1:60){
    caRNA_Slope_list<-rep(0,2418*4)
    for(i in 1:2418){
      ind_delta<-c((15-slope_delta),(15+slope_delta),(30-slope_delta),(30+slope_delta),(60-slope_delta),(60+slope_delta))
      ind_delta[ind_delta<0]<-0
      ind_delta[ind_delta>80]<-80
      temp<-pchip(c(0,15,30,60),caRNA_raw.mean[(i*4-3):(i*4)],ind_delta)
      
      caRNA_Slope_list[(i*4-2):(i*4)]<-c((temp[2]-temp[1]),(temp[4]-temp[3]),(temp[6]-temp[5]))/(c((ind_delta[2]-ind_delta[1]),(ind_delta[4]-ind_delta[3]),(ind_delta[6]-ind_delta[5])))
      #slope for non activated to be 0
      if(all(caRNA_raw.mean[(i*4-3):(i*4)]<3)){
        caRNA_Slope_list[(i*4-2):(i*4)]<-0
      }
    }
    #set negative slope to be 0
    caRNA_Slope_list[caRNA_Slope_list<0]<-0
    caRNA_Slope_Full_list<-rbind(caRNA_Slope_Full_list,caRNA_Slope_list)
  }
  return(caRNA_Slope_Full_list)
}
########################################################################
# estimated_var<-Prior_TV_bio_V0.01_BT50_rep1000
# ground_truth_var<-var_empirical_bio_V0.01_BT50_rep1000
# plot_name<-"Prior_TV_estimate_pure_bio_V0.01_BT50_rep1000"
plot_heatmap<-function(ground_truth_var,estimated_var,plot_name){
  x1<-matrix(t(ground_truth_var),ncol=1)
  x2<-estimated_var$var_est_each
  df <- data.frame(x1,x2)
  df$density <- get_density(df$x1, df$x2)
  
  cor_cal<-cor(df$x1, df$x2)
  my_text <- paste("Pearson Correlation",round(cor_cal,digits = 2))
  my_grob = grid.text(my_text, x=0.45,  y=0.85, gp=gpar(col="firebrick", fontsize=15, fontface="bold"))
  
  ggplot(df) + geom_point(aes(x1, x2, color = density),size = 1) + scale_color_viridis()+
    labs(x = "Ground Truth Variance",y = "Estimated Variance")+
    xlim(0, 2000)+ 
    ylim(0, 2000)+ 
    theme_bw()+
    annotation_custom(my_grob)+
    geom_abline(color='red',slope=1,size=1.5,linetype = "dashed")+
    theme(axis.text.x=element_text(size=rel(1.8),face="bold"),axis.text.y=element_text(size=rel(1.5),face="bold"),
          axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
          legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
    ggsave(paste(plot_name,"_density.pdf"),width=8,height=6)
}
plot_heatmap_N<-function(ground_truth_var,estimated_var,plot_name){
  x1<-matrix(t(ground_truth_var),ncol=1)
  x2<-estimated_var$var_est_each
  df <- data.frame(x1,x2)
  df$density <- get_density(df$x1, df$x2)
  
  cor_cal<-cor(df$x1, df$x2)
  my_text <- paste("Pearson Correlation",round(cor_cal,digits = 2))
  my_grob = grid.text(my_text, x=0.45,  y=0.85, gp=gpar(col="firebrick", fontsize=15, fontface="bold"))
  
  ggplot(df) + geom_point(aes(x1, x2, color = density),size = 1) + scale_color_viridis()+
    labs(x = " ",y = " ")+
    xlim(0, 3000)+ 
    ylim(0, 3000)+ 
    theme_bw()+
    annotation_custom(my_grob)+ 
    geom_abline(color='red',slope=1,size=1.5,linetype = "dashed")+
    theme(axis.text.x=element_text(size=rel(1.8),face="bold"),axis.text.y=element_text(size=rel(1.5),face="bold"),
          axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
          legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
    ggsave(paste(plot_name,"_density.pdf"),width=8,height=6)
}
#############################Load Data V4#######################################
################################################################################
load("Nondyn_active_temp_BT50_parV_rep1000_N_M0.75.RData")
load("Nondyn_bio_V0.01_BT50_rep1000_M0.75.RData")
load("var_empirical_pure_bio_V0.01_BT50_rep1000_M0.75.RData")
load("var_empirical_bio_V0.01_BT50_rep1000_M0.75.RData")

load("Nondyn_active_GT.RData")
load("slope_list_GT.RData")

caRNA_raw_BT50_parV_rep1000<-cbind(matrix(t(Nondyn_active_temp_BT50_parV_rep1000_N_M0.75[1:93,]),ncol=1),matrix(t(Nondyn_active_temp_BT50_parV_rep1000_N_M0.75[94:186,]),ncol=1))
caRNA_raw_BT50_parV_rep1000.mean<-apply(caRNA_raw_BT50_parV_rep1000,1,mean)

caRNA_raw_bio_V0.01_BT50_rep1000<-cbind(matrix(t(Nondyn_bio_V0.01_BT50_rep1000_M0.75[,1:104]),ncol=1),matrix(t(Nondyn_bio_V0.01_BT50_rep1000_M0.75[,105:208]),ncol=1))
caRNA_raw_bio_V0.01_BT50_rep1000.mean<-apply(caRNA_raw_bio_V0.01_BT50_rep1000,1,mean)
####################################################################
###########Slope Estimate and Plot
####################################################################
caRNA_Slope_BT50_parV_rep1000<-Slope_estimate_f(caRNA_raw_BT50_parV_rep1000.mean)
caRNA_Slope_bio_V0.01_BT50_rep1000<-Slope_estimate_f(caRNA_raw_bio_V0.01_BT50_rep1000.mean)
save(caRNA_Slope_BT50_parV_rep1000,file="caRNA_Slope_BT50_parV_rep1000.RData")
save(caRNA_Slope_bio_V0.01_BT50_rep1000,file="caRNA_Slope_bio_V0.01_BT50_rep1000.RData")
# load("caRNA_Slope_BT50_parV_rep1000.RData")
# load("caRNA_Slope_bio_V0.01_BT50_rep1000.RData")
################ggplot density
################only para vary
slope_ind<-round(2*0.75*(50)^0.5)
if(slope_ind>60) slope_ind<-60
if(slope_ind<1) slope_ind<-1

x1<-slope_list_GT[slope_ind,]
x2<-caRNA_Slope_BT50_parV_rep1000[slope_ind,]
df_MAP <- data.frame(x1,x2)
df_MAP$density <- get_density(df_MAP$x1, df_MAP$x2)

Slope_cor<-cor(df_MAP$x1, df_MAP$x2)
my_text <- paste("Pearson Correlation",round(Slope_cor,digits = 2))
my_grob = grid.text(my_text, x=0.25,  y=0.85, gp=gpar(col="firebrick", fontsize=15, fontface="bold"))

ggplot(df_MAP) + geom_point(aes(x1, x2, color = density),size = 1) + scale_color_viridis()+
  labs(x = "Ground Truth Slope",y = "Estimated Slope")+
  theme_bw()+
  annotation_custom(my_grob)+
  geom_abline(color='red',slope=1,size=1.5,linetype = "dashed")+
  theme(axis.text.x=element_text(size=rel(1.8),face="bold"),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Slope_estimation_BT50_parV_rep1000_density_M0.75up.pdf",width=8,height=6)
################all noise
x1<-slope_list_GT[slope_ind,]
x2<-caRNA_Slope_bio_V0.01_BT50_rep1000[slope_ind,]
df_MAP <- data.frame(x1,x2)
df_MAP$density <- get_density(df_MAP$x1, df_MAP$x2)

Slope_cor<-cor(df_MAP$x1, df_MAP$x2)
my_text <- paste("Pearson Correlation",round(Slope_cor,digits = 2))
my_grob = grid.text(my_text, x=0.25,  y=0.85, gp=gpar(col="firebrick", fontsize=15, fontface="bold"))

ggplot(df_MAP) + geom_point(aes(x1, x2, color = density),size = 1) + scale_color_viridis()+
  labs(x = "Ground Truth Slope",y = "Estimated Slope")+
  theme_bw()+
  annotation_custom(my_grob)+
  geom_abline(color='red',slope=1,size=1.5,linetype = "dashed")+
  theme(axis.text.x=element_text(size=rel(1.8),face="bold"),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Slope_estimation_bio_V0.01_BT50_rep1000_density_M0.75up.pdf",width=8,height=6)
####################################################################
#######################################Variance Estimation############################
#####################################all parameters###############################
Prior_Con_bio_V0.01_BT50_rep1000_M0.75up<-Prior_Con_estimate(caRNA_raw_bio_V0.01_BT50_rep1000)
save(Prior_Con_bio_V0.01_BT50_rep1000_M0.75up,file = "Prior_Con_bio_V0.01_BT50_rep1000_M0.75up.RData")

Prior_TV_bio_V0.01_BT50_rep1000_M0.75up<-Prior_TV_estimate(caRNA_raw_bio_V0.01_BT50_rep1000,caRNA_Slope_bio_V0.01_BT50_rep1000)
save(Prior_TV_bio_V0.01_BT50_rep1000_M0.75up,file = "Prior_TV_bio_V0.01_BT50_rep1000_M0.75up.RData")
########################################
plot_heatmap(ground_truth_var=var_empirical_bio_V0.01_BT50_rep1000_M0.75,estimated_var=Prior_TV_bio_V0.01_BT50_rep1000_M0.75up,
             plot_name="Prior_TV_estimate_all_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_bio_V0.01_BT50_rep1000_M0.75,estimated_var=Prior_Con_bio_V0.01_BT50_rep1000_M0.75up,
             plot_name="Prior_Con_estimate_all_bio_V0.01_BT50_rep1000_M0.75up" )
########################################
MAP_TV_bio_V0.01_BT50_rep1000_M0.75up<-MAP_TV_estimate(caRNA_raw_bio_V0.01_BT50_rep1000,caRNA_Slope_bio_V0.01_BT50_rep1000,var_trend_n=Prior_TV_bio_V0.01_BT50_rep1000_M0.75up$var_est_n,
                                                       var_trend_g=Prior_TV_bio_V0.01_BT50_rep1000_M0.75up$var_est_g,0.01)
save(MAP_TV_bio_V0.01_BT50_rep1000_M0.75up,file = "MAP_TV_bio_V0.01_BT50_rep1000_M0.75up.RData")

MAP_Con_bio_V0.01_BT50_rep1000_M0.75up<-MAP_Con_estimate(caRNA_raw_bio_V0.01_BT50_rep1000,var_trend_n=Prior_Con_bio_V0.01_BT50_rep1000_M0.75up$var_est_n,
                                                         var_trend_g=Prior_Con_bio_V0.01_BT50_rep1000_M0.75up$var_est_g,0.01)
save(MAP_Con_bio_V0.01_BT50_rep1000_M0.75up,file = "MAP_Con_bio_V0.01_BT50_rep1000_M0.75up.RData")

MLE_TV_bio_V0.01_BT50_rep1000_M0.75up<-MLE_TV_estimate(caRNA_raw_bio_V0.01_BT50_rep1000,caRNA_Slope_bio_V0.01_BT50_rep1000)
save(MLE_TV_bio_V0.01_BT50_rep1000_M0.75up,file = "MLE_TV_bio_V0.01_BT50_rep1000_M0.75up.RData")
MLE_Con_bio_V0.01_BT50_rep1000_M0.75up<-MLE_Con_estimate(caRNA_raw_bio_V0.01_BT50_rep1000)
save(MLE_Con_bio_V0.01_BT50_rep1000_M0.75up,file = "MLE_Con_bio_V0.01_BT50_rep1000_M0.75up.RData")
################################################################################################################
########################################pure parameter vary########################################################
################################################################################################################
Prior_TV_BT50_parV_rep1000_M0.75up<-Prior_TV_estimate(caRNA_raw_BT50_parV_rep1000,caRNA_Slope_BT50_parV_rep1000)
save(Prior_TV_BT50_parV_rep1000_M0.75up,file = "Prior_TV_BT50_parV_rep1000_M0.75up.RData")
Prior_Con_BT50_parV_rep1000_M0.75up<-Prior_Con_estimate(caRNA_raw_BT50_parV_rep1000)
save(Prior_Con_BT50_parV_rep1000_M0.75up,file = "Prior_Con_BT50_parV_rep1000_M0.75up.RData")

MLE_TV_BT50_parV_rep1000_M0.75up<-MLE_TV_estimate(caRNA_raw_BT50_parV_rep1000,caRNA_Slope_BT50_parV_rep1000)
save(MLE_TV_BT50_parV_rep1000_M0.75up,file = "MLE_TV_BT50_parV_rep1000_M0.75up.RData")
MLE_Con_BT50_parV_rep1000_M0.75up<-MLE_Con_estimate(caRNA_raw_BT50_parV_rep1000)
save(MLE_Con_BT50_parV_rep1000_M0.75up,file = "MLE_Con_BT50_parV_rep1000_M0.75up.RData")

MAP_TV_BT50_parV_rep1000_M0.75up<-MAP_TV_estimate(caRNA_raw_BT50_parV_rep1000,caRNA_Slope_BT50_parV_rep1000,var_trend_n=Prior_TV_BT50_parV_rep1000_M0.75up$var_est_n,
                                                  var_trend_g=Prior_TV_BT50_parV_rep1000_M0.75up$var_est_g, 0.01)
save(MAP_TV_BT50_parV_rep1000_M0.75up,file = "MAP_TV_BT50_parV_rep1000_M0.75up.RData")
MAP_Con_BT50_parV_rep1000_M0.75up<-MAP_Con_estimate(caRNA_raw_BT50_parV_rep1000,var_trend_n=Prior_Con_BT50_parV_rep1000_M0.75up$var_est_n,
                                                    var_trend_g=Prior_Con_BT50_parV_rep1000_M0.75up$var_est_g, 0.01)
save(MAP_Con_BT50_parV_rep1000_M0.75up,file = "MAP_Con_BT50_parV_rep1000_M0.75up.RData")
#############################################################################################
#############################################################################################
# load("MLE_TV_BT50_parV_rep1000_M0.75up.RData")
# load("MLE_Con_BT50_parV_rep1000_M0.75up.RData")
# 
# load("Prior_TV_BT50_parV_rep1000_M0.75up.RData")
# load("Prior_Con_BT50_parV_rep1000_M0.75up.RData")
# 
# load("MAP_TV_BT50_parV_rep1000_M0.75up.RData")
# load("MAP_Con_BT50_parV_rep1000_M0.75up.RData")
# 
# load("MLE_TV_bio_V0.01_BT50_rep1000_M0.75up.RData")
# load("MLE_Con_bio_V0.01_BT50_rep1000_M0.75up.RData")
# 
# load("Prior_TV_bio_V0.01_BT50_rep1000_M0.75up.RData")
# load("Prior_Con_bio_V0.01_BT50_rep1000_M0.75up.RData")
# 
# load("MAP_TV_bio_V0.01_BT50_rep1000_M0.75up.RData")
# load("MAP_Con_bio_V0.01_BT50_rep1000_M0.75up.RData")
#################################
#################Plot#################
##################################
###########################Para only Biology
plot_heatmap(ground_truth_var=var_empirical_pure_bio_V0.01_BT50_rep1000_M0.75,estimated_var=MLE_TV_BT50_parV_rep1000_M0.75up,
             plot_name="MLE_TV_estimate_pure_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_pure_bio_V0.01_BT50_rep1000_M0.75,estimated_var=MLE_Con_BT50_parV_rep1000_M0.75up,
             plot_name="MLE_Con_estimate_pure_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_pure_bio_V0.01_BT50_rep1000_M0.75,estimated_var=Prior_TV_BT50_parV_rep1000_M0.75up,
             plot_name="Prior_TV_estimate_pure_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_pure_bio_V0.01_BT50_rep1000_M0.75,estimated_var=Prior_Con_BT50_parV_rep1000_M0.75up,
             plot_name="Prior_Con_estimate_pure_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_pure_bio_V0.01_BT50_rep1000_M0.75,estimated_var=MAP_TV_BT50_parV_rep1000_M0.75up,
             plot_name="MAP_TV_estimate_pure_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_pure_bio_V0.01_BT50_rep1000_M0.75,estimated_var=MAP_Con_BT50_parV_rep1000_M0.75up,
             plot_name="MAP_Con_estimate_pure_bio_V0.01_BT50_rep1000_M0.75up" )
###########################All Biology
plot_heatmap(ground_truth_var=var_empirical_bio_V0.01_BT50_rep1000_M0.75,estimated_var=MLE_TV_bio_V0.01_BT50_rep1000_M0.75up,
             plot_name="MLE_TV_estimate_all_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_bio_V0.01_BT50_rep1000_M0.75,estimated_var=MLE_Con_bio_V0.01_BT50_rep1000_M0.75up,
             plot_name="MLE_Con_estimate_all_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_bio_V0.01_BT50_rep1000_M0.75,estimated_var=Prior_TV_bio_V0.01_BT50_rep1000_M0.75up,
             plot_name="Prior_TV_estimate_all_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_bio_V0.01_BT50_rep1000_M0.75,estimated_var=Prior_Con_bio_V0.01_BT50_rep1000_M0.75up,
             plot_name="Prior_Con_estimate_all_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_bio_V0.01_BT50_rep1000_M0.75,estimated_var=MAP_TV_bio_V0.01_BT50_rep1000_M0.75up,
             plot_name="MAP_TV_estimate_all_bio_V0.01_BT50_rep1000_M0.75up" )

plot_heatmap(ground_truth_var=var_empirical_bio_V0.01_BT50_rep1000_M0.75,estimated_var=MAP_Con_bio_V0.01_BT50_rep1000_M0.75up,
             plot_name="MAP_Con_estimate_all_bio_V0.01_BT50_rep1000_M0.75up" )
