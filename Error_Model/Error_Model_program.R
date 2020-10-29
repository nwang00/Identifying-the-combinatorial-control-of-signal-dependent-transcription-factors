rm(list=ls())
library(signal)
library(MASS)
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
Prior_TV_estimate<-function(data.raw,caRNA_Slope_Full_list,Rep_num){
  # data.raw<-Gene_exp_raw
  # caRNA_Slope_Full_list<-Gene_exp.slope
  data.mean<-apply(data.raw,1,mean)
  
  normal_ind<-which(data.mean>3)
  zero_ind<-which(apply(data.raw,1,min)==0)
  gamma_ind<-setdiff(setdiff(c(1:length(data.mean)),normal_ind),zero_ind)
  
  gamma_ind_all<-rep(gamma_ind, Rep_num)+rep(c(0:(Rep_num-1))*dim(data.raw)[1],each=length(gamma_ind))
  normal_ind_all<-rep(normal_ind, Rep_num)+rep(c(0:(Rep_num-1))*dim(data.raw)[1],each=length(normal_ind))

  K1p<-runif(1000, min=-20, max=10)
  K2p<-runif(1000, min=-20, max=10)
  K3p<-runif(1000, min=-20, max=10)
  K4p<-runif(1000, min=-20, max=10)
  
  RNA_gw_list<-vector()
  print("Start estimation normal")
  pb = txtProgressBar(min = 1, max = 200)
  
  for(i in 1:200){
    RNA_LH<-optim(c(K1p[i],K2p[i],K3p[i],K4p[i]),method="L-BFGS-B",fn=NLL_Prior_TV_n,lower=rep((-20),4),upper=rep((10),4),normal_ind_all=normal_ind_all,
                  data_point=matrix(data.raw,ncol=1),slop_square=((caRNA_Slope_Full_list)^2)[,rep(c(1:length(data.mean)),Rep_num)],
                  mean_est=rep(data.mean,Rep_num))
    
    RNA_gw_list<-c(RNA_gw_list,RNA_LH)
    setTxtProgressBar(pb,i)
  }
  
  var_est_n<-2*10^(RNA_gw_list[order(unlist(RNA_gw_list[seq(2,5*1000,5)]))[1]*5-4]$par)
  ###############################################################
  ###############################################################
  K1p<-runif(1000, min=-20, max=10)
  
  RNA_gw_list<-vector()
  print("Start estimation gamma")
  pb = txtProgressBar(min = 1, max = 200)
  
  for(i in 1:200){
    RNA_LH<-optim(c(K1p[i]),method="L-BFGS-B",fn=NLL_Prior_g,lower=(-20),upper=(10),gamma_ind_all=gamma_ind_all,
                  data_point=matrix(data.raw,ncol=1),
                  mean_est=rep(data.mean,Rep_num))
    
    RNA_gw_list<-c(RNA_gw_list,RNA_LH)
    setTxtProgressBar(pb,i)
  }
  
  var_est_g<-2*10^(RNA_gw_list[order(unlist(RNA_gw_list[seq(2,5*1000,5)]))[1]*5-4]$par)
  ###############################################################
  return(list(data.mean=data.mean,var_est_g=var_est_g,var_est_n=var_est_n))
}
#################################Slope Estimate######################
# caRNA_raw.mean<-Gene_exp.mean
# Time_course_slp<-Time_course
Slope_estimate_f<-function(caRNA_raw.mean,Time_course_slp){
  max_range<-max(Time_course_slp)
  caRNA_Slope_Full_list<-vector()
  for(slope_delta in 1:max_range){
    
    caRNA_Slope_list<-rep(0,length(caRNA_raw.mean))
    
    for(i in 1:(length(caRNA_raw.mean)/length(Time_course_slp))){
      
      ind_delta<-rep(Time_course_slp,each=2)+rep(c(-slope_delta,slope_delta),length(Time_course_slp))

      ind_delta[ind_delta<0]<-0
      ind_delta[ind_delta>(max_range+20)]<-(max_range+20)
      
      temp<-pchip(Time_course_slp,caRNA_raw.mean[((i-1)*length(Time_course_slp)+1):(i*length(Time_course_slp))],ind_delta)
      
      caRNA_Slope_list[((i-1)*length(Time_course_slp)+1):(i*length(Time_course_slp))]<-(temp[seq(2,(length(Time_course_slp)*2),2)]-temp[seq(1,(length(Time_course_slp)*2),2)])/(ind_delta[seq(2,(length(Time_course_slp)*2),2)]-ind_delta[seq(1,(length(Time_course_slp)*2),2)])

      caRNA_Slope_list[((i-1)*length(Time_course_slp)+which(Time_course_slp==0))]<-0
      
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
#############################Read Parameters#######################################
################################################################################
Para_file_path<- commandArgs(trailingOnly=TRUE)
Paras_file<-read.table(file=Para_file_path)
# Paras_file<-read.table(file="./Parameters_ErrorModel.txt")
Time_course<-as.numeric(unlist(strsplit(Paras_file[1,2],",")))
Rep_num<-as.numeric(Paras_file[2,2])
Gene_exp_matrix<-read.table(file=Paras_file[3,2],header = T)
########################################################################
#############################Calculate Mean#######################################
################################################################################
Gene_exp_raw<-vector()
for(i in 1:Rep_num){
  Gene_exp_raw<-cbind(Gene_exp_raw,matrix(t(Gene_exp_matrix[,seq((i+1),dim(Gene_exp_matrix)[2],Rep_num)]),ncol=1))
}

Gene_exp.mean<-apply(Gene_exp_raw,1,mean)
########################################################################
#############################Calculate Slope#######################################
################################################################################
Gene_exp.slope<-Slope_estimate_f(Gene_exp.mean,Time_course)
########################################################################
#############################Calculate Variance#######################################
################################################################################
GeneExp.Variance<-Prior_TV_estimate(Gene_exp_raw,Gene_exp.slope,Rep_num)

print(paste("Basal_variance",round(GeneExp.Variance$var_est_g,digits = 5)))
print(paste("Inducible_variance_parameters(x0)",round(GeneExp.Variance$var_est_n[1],digits = 5)))
print(paste("Inducible_variance_parameters(x1)",round(GeneExp.Variance$var_est_n[2],digits = 5)))
print(paste("Inducible_variance_parameters(x2)",round(GeneExp.Variance$var_est_n[3],digits = 5)))
print(paste("Temporal_varaince",round(GeneExp.Variance$var_est_n[4],digits = 5)))

write.table(c(GeneExp.Variance$var_est_g,GeneExp.Variance$var_est_n),file = "Inferred_uncertainly_level.txt",row.names = c("Basal_variance","Inducible_variance_parameters(x0)","Inducible_variance_parameters(x1)","Inducible_variance_parameters(x2)","Temporal_varaince"),
            quote = F,col.names = F)
save(GeneExp.Variance,file = "GeneExp.Variance.RData")
