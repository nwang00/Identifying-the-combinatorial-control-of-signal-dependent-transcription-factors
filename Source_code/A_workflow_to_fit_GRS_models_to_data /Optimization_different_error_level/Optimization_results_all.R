rm(list=ls())
library(ggplot2)
library(scales)
library(prodlim)
source("Nondyn_opt_newbasalnew.R")
###########################################################################################
load("Nondyn_WT_bio.RData")
Nonredundent_name<-read.table(file = "Nonredundent_name.txt",sep = "\t")
logic_name<-c("1a2a3","1o2o3","1o(2a3)","1a(2o3)","2o(1a3)","2a(1o3)","3o(1a2)","3a(1o2)")
para_name<-paste(expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,1],
                 expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,2],
                 expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,3],sep="")[-27]
GRS_name<-paste(rep(logic_name,each=26),rep(para_name,8))
para_list<-as.data.frame(expand.grid(c(0.1,1,10),c(0.1,1,10),c(0.1,1,10)))[-27,]
para_all<-cbind(para_list[rep(c(1:26),8),],1,0.1,0)
para_true<-para_all[which(GRS_name%in%apply(Nonredundent_name,1,as.character)),]
para_true[,4]<-Nondyn_WT_bio
###########################################################################################
# load("Nondyn_raw_M1_2.RData")
# load("Nondyn_bio_V0.01_B5T5_rep1000_M1_2.RData")

# data_name<-"Nondyn_bio_V0.01_B5T5_rep1000_0.02"
# simulate_name<-"Nondyn_raw_0.02"

para_estimator<-function(data_name,simulate_name){
  
  load(paste(data_name,".RData",sep=""))
  load(paste(simulate_name,".RData",sep=""))
  ###########################################################################################
  simulate_data<-get(data_name)
  ###########################################################################################
  mse_GRNs<-vector()
  par_GRNs<-vector()
  for(args in 1:93){
    run_id<-rownames(simulate_data)[as.numeric(args)]
    name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
    for(GRN_ind in 1:8){
      sim_res_GRS<-get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))
      min_ind<-which(unlist(sim_res_GRS[seq(2,1800,6)])==min(unlist(sim_res_GRS[seq(2,1800,6)])))
      mse_GRNs<-rbind(mse_GRNs,unlist(sim_res_GRS[min_ind*6-4]))
      par_GRNs[length(par_GRNs)+1]<-list(unlist(get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))[min_ind*6-5]))
    }
    print(args)
  }
  sim_GRS<-t(matrix(mse_GRNs,nrow=8))
  ####Identify select model
  sim_GRS_min<-matrix(0,ncol = 8,nrow=93)
  for(i in 1:93){
    sim_GRS_min[i,which(sim_GRS[i,]==min(sim_GRS[i,]))]<-1
  }
  ####Estimated GRS
  min_list<-apply(sim_GRS,1,function(x) which(x==min(x)))
  ####ground truth index
  min_list_GT<-c(rep(1,7),rep(2,26),rep(3,12),rep(4,8),rep(5,12),rep(6,8),rep(7,12),rep(8,8))
  ####Likelihood ratio
  sim_GRS_F<-rep(0,93)
  for(i in 1:93){
    sim_GRS_F[i]<-log2(exp(1))*(min(sim_GRS[i,-min_list_GT[i]])-sim_GRS[i,min_list_GT[i]])
  }
  ###########################################################################################
  para_est<-t(matrix(unlist(par_GRNs[min_list+seq(0,92*8,8)]),nrow=6))
  para_fold<-vector()
  for(i in 1:93){
    para_fold<-rbind(para_fold,c(10^para_est[i,1:5])/para_true[i,1:5])
  } 
  para_fold<-cbind(para_fold,para_est[,6])
  
  return(para_fold)
}

para_fold_raw_0.02<-para_estimator("Nondyn_bio_V0.01_B5T5_rep1000_0.02","Nondyn_raw_0.02")
para_fold_raw_M1_2<-para_estimator("Nondyn_bio_V0.01_B5T5_rep1000_M1_2","Nondyn_raw_M1_2")
para_fold_raw_M0.75<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M0.75","Nondyn_raw_M0.75")
para_fold_raw_M1.05<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M1.05","Nondyn_raw_M1.05")
para_fold_raw_M1.2<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M1.2","Nondyn_raw_M1.2")

para_fold_TV_0.02<-para_estimator("Nondyn_bio_V0.01_B5T5_rep1000_0.02","Nondyn_V0.01_TB5_TT5_0.02_fupdate")
para_fold_TV_M1_2<-para_estimator("Nondyn_bio_V0.01_B5T5_rep1000_M1_2","Nondyn_V0.01_BT50_M1_2_fupdate")
para_fold_TV_M0.75<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M0.75","Nondyn_V0.01_BT50_M0.75_fupdate")
para_fold_TV_M1.05<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M1.05","Nondyn_V0.01_BT50_M1.05_fupdate")
para_fold_TV_M1.2<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M1.2","Nondyn_V0.01_BT50_M1.2_fupdate")

para_fold_Con_0.02<-para_estimator("Nondyn_bio_V0.01_B5T5_rep1000_0.02","Nondyn_Con_V0.01_TB5_TT5_0.02_fupdate")
para_fold_Con_M1_2<-para_estimator("Nondyn_bio_V0.01_B5T5_rep1000_M1_2","Nondyn_Con_V0.01_BT50_M1_2_fupdate")
para_fold_Con_M0.75<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M0.75","Nondyn_Con_V0.01_BT50_M0.75_fupdate")
para_fold_Con_M1.05<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M1.05","Nondyn_Con_V0.01_BT50_M1.05_fupdate")
para_fold_Con_M1.2<-para_estimator("Nondyn_bio_V0.01_BT50_rep1000_M1.2","Nondyn_Con_V0.01_BT50_M1.2_fupdate")

plot_Kds<-function(para_fold_TV,para_fold_Con,para_fold_raw,plot_name){
  para_est_plot_Kds<-data.frame(para=rep(c("Kd1_TV","Kd2_TV","Kd3_TV","Kd1_Con","Kd2_Con","Kd3_Con","Kd1_raw","Kd2_raw","Kd3_raw"),each=93),est_rel=c(para_fold_TV[,1],para_fold_TV[,2],para_fold_TV[,3],para_fold_Con[,1],para_fold_Con[,2],para_fold_Con[,3],para_fold_raw[,1],para_fold_raw[,2],para_fold_raw[,3]))
  para_est_plot_Kds$para<-factor(para_est_plot_Kds$para, levels = c("Kd1_TV","Kd1_Con","Kd1_raw","Kd2_TV","Kd2_Con","Kd2_raw","Kd3_TV","Kd3_Con","Kd3_raw"))
  #########################################
  pp <- ggplot(para_est_plot_Kds, aes(x=para,y=est_rel)) + geom_boxplot(mapping=aes(fill=para),outlier.shape = NA,color=rep(c("red","black","deepskyblue3"),3),alpha = 0.3)
  pp+labs(x="",y = "")+
    geom_hline(yintercept=1,color = adjustcolor("orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
    scale_y_continuous(trans='log10',limits = c(0.01,100))+
    scale_fill_manual(values=rep(c("red","black","deepskyblue3"),3))+
    geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),aes(color = para))+
    scale_color_manual(values=rep(c("red","black","deepskyblue3"),3))+
    scale_x_discrete(labels=rep(c("Kd1","Kd2","Kd3"),each=3))+
    theme_bw()+
    theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
          axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
          legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
    ggsave(paste(plot_name,".pdf",sep=""),width=6,height=6)
}

plot_Kds(para_fold_TV_0.02,para_fold_Con_0.02,para_fold_raw_0.02,"Para_est_Kds_0.02_fupdate")
plot_Kds(para_fold_TV_M1_2,para_fold_Con_M1_2,para_fold_raw_M1_2,"Para_est_Kds_M1_2_fupdate")
plot_Kds(para_fold_TV_M0.75,para_fold_Con_M0.75,para_fold_raw_M0.75,"Para_est_Kds_M0.75_fupdate")
plot_Kds(para_fold_TV_M1.05,para_fold_Con_M1.05,para_fold_raw_M1.05,"Para_est_Kds_M1.05_fupdate")
plot_Kds(para_fold_TV_M1.2,para_fold_Con_M1.2,para_fold_raw_M1.2,"Para_est_Kds_M1.2_fupdate")
