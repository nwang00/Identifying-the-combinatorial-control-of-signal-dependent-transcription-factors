setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics")
rm(list=ls())
library(ggplot2)
library(scales)
library(prodlim)

load("Nondyn_bio_V0.01_B5T5_rep1000_0.02.RData")
load("Nondyn_V0.01_TB5_TT5_0.02_fupdate.RData")
source("Nondyn_opt_newbasalnew.R")

mse_GRNs_TV<-vector()
par_GRNs_TV<-vector()
for(args in 1:93){
  run_id<-rownames(Nondyn_bio_V0.01_B5T5_rep1000_0.02)[as.numeric(args)]
  name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
  for(GRN_ind in 1:8){
    sim_res_GRS<-get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))
    min_ind<-which(unlist(sim_res_GRS[seq(2,1800,6)])==min(unlist(sim_res_GRS[seq(2,1800,6)])))
    mse_GRNs_TV<-rbind(mse_GRNs_TV,unlist(sim_res_GRS[min_ind*6-4]))
    par_GRNs_TV[length(par_GRNs_TV)+1]<-list(unlist(get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))[min_ind*6-5]))
  }
  print(args)
}
sim_GRS_TV<-t(matrix(mse_GRNs_TV,nrow=8))
####Identify select model
sim_GRS_min_TV<-matrix(0,ncol = 8,nrow=93)
for(i in 1:93){
  sim_GRS_min_TV[i,which(sim_GRS_TV[i,]==min(sim_GRS_TV[i,]))]<-1
}
####Estimated GRS
min_list_TV<-apply(sim_GRS_TV,1,function(x) which(x==min(x)))
####ground truth index
min_list_GT<-c(rep(1,7),rep(2,26),rep(3,12),rep(4,8),rep(5,12),rep(6,8),rep(7,12),rep(8,8))
####Likelihood ratio
sim_GRS_F_TV<-rep(0,93)
for(i in 1:93){
  sim_GRS_F_TV[i]<-log2(exp(1))*(min(sim_GRS_TV[i,-min_list_GT[i]])-sim_GRS_TV[i,min_list_GT[i]])
}
#################################################################
################Converge Test####################################
#################################################################
num_list<-seq(100,300,20)
GRS_F_list<-vector()
for(args in 1:93){
  run_id<-rownames(Nondyn_bio_V0.01_B5T5_rep1000_0.02)[as.numeric(args)]
  name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
  
  GRS_F_itm<-vector()
  for(i in 1:length(num_list)){
    mse_GRNs<-vector()
    for(GRN_ind in 1:8){
      sim_res_GRS<-get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))
      min_ind<-which(unlist(sim_res_GRS[seq(2,6*num_list[i],6)])==min(unlist(sim_res_GRS[seq(2,6*num_list[i],6)])))
      mse_GRNs<-c(mse_GRNs,unlist(sim_res_GRS[min_ind*6-4]))
    }
    GRS_F_itm<-cbind(GRS_F_itm,log2(exp(1))*(min(mse_GRNs[-min_list_GT[args]])-mse_GRNs[min_list_GT[args]]))
  }
  
  GRS_F_list<-rbind(GRS_F_list,GRS_F_itm)
}

Converge_plot<-data.frame(Ite_num=rep(num_list,each=93),Converge=t(matrix(GRS_F_list,nrow=1)))
colnames(Converge_plot)[1:2]<-c("Ite_num","Converge")
Converge_plot$Ite_num<-factor(Converge_plot$Ite_num, levels = num_list)

# pp <- ggplot(Converge_plot, aes(x=Ite_num,y=Converge)) + geom_boxplot(mapping=aes(fill=Ite_num),color=rep("black",length(num_list)),outlier.shape = NA,alpha = 0.3)
# pp+labs(x="",y = "")+
#   scale_fill_manual(values=rep("black",length(num_list)))+
#   geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2))+
#   scale_color_manual(values=rep("black",length(num_list)))+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
#         axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
#         legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
#   ggsave("Converge_TV_0.02_fupdate.pdf",width=10,height=6)
###########################################################################################
#############################Conventional##################################################
###########################################################################################
load("Nondyn_Con_V0.01_TB5_TT5_0.02_fupdate.RData")

mse_GRNs_Con<-vector()
par_GRNs_Con<-vector()
for(args in 1:93){
  run_id<-rownames(Nondyn_bio_V0.01_B5T5_rep1000_0.02)[as.numeric(args)]
  name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
  for(GRN_ind in 1:8){
    sim_res_GRS<-get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))
    min_ind<-which(unlist(sim_res_GRS[seq(2,1800,6)])==min(unlist(sim_res_GRS[seq(2,1800,6)])))
    mse_GRNs_Con<-rbind(mse_GRNs_Con,unlist(sim_res_GRS[min_ind*6-4]))
    par_GRNs_Con[length(par_GRNs_Con)+1]<-list(unlist(get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))[min_ind*6-5]))
  }
  print(args)
}
sim_GRS_Con<-t(matrix(mse_GRNs_Con,nrow=8))
####Identify select model
sim_GRS_min_Con<-matrix(0,ncol = 8,nrow=93)
for(i in 1:93){
  sim_GRS_min_Con[i,which(sim_GRS_Con[i,]==min(sim_GRS_Con[i,]))]<-1
}
####Estimated GRS
min_list_Con<-apply(sim_GRS_Con,1,function(x) which(x==min(x)))
####ground truth index
min_list_GT<-c(rep(1,7),rep(2,26),rep(3,12),rep(4,8),rep(5,12),rep(6,8),rep(7,12),rep(8,8))
####Likelihood ratio
sim_GRS_F_Con<-rep(0,93)
for(i in 1:93){
  sim_GRS_F_Con[i]<-log2(exp(1))*(min(sim_GRS_Con[i,-min_list_GT[i]])-sim_GRS_Con[i,min_list_GT[i]])
}
#################################################################
################Converge Test####################################
#################################################################
num_list<-seq(100,300,20)
GRS_F_list<-vector()
for(args in 1:93){
  run_id<-rownames(Nondyn_bio_V0.01_B5T5_rep1000_0.02)[as.numeric(args)]
  name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
  
  GRS_F_itm<-vector()
  for(i in 1:length(num_list)){
    mse_GRNs<-vector()
    for(GRN_ind in 1:8){
      sim_res_GRS<-get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))
      min_ind<-which(unlist(sim_res_GRS[seq(2,6*num_list[i],6)])==min(unlist(sim_res_GRS[seq(2,6*num_list[i],6)])))
      mse_GRNs<-c(mse_GRNs,unlist(sim_res_GRS[min_ind*6-4]))
    }
    GRS_F_itm<-cbind(GRS_F_itm,log2(exp(1))*(min(mse_GRNs[-min_list_GT[args]])-mse_GRNs[min_list_GT[args]]))
  }
  
  GRS_F_list<-rbind(GRS_F_list,GRS_F_itm)
}

Converge_plot<-data.frame(Ite_num=rep(num_list,each=93),Converge=t(matrix(GRS_F_list,nrow=1)))
colnames(Converge_plot)[1:2]<-c("Ite_num","Converge")
Converge_plot$Ite_num<-factor(Converge_plot$Ite_num, levels = num_list)

# pp <- ggplot(Converge_plot, aes(x=Ite_num,y=Converge)) + geom_boxplot(mapping=aes(fill=Ite_num),color=rep("black",length(num_list)),outlier.shape = NA,alpha = 0.3)
# pp+labs(x="",y = "")+
#   scale_fill_manual(values=rep("black",length(num_list)))+
#   geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2))+
#   scale_color_manual(values=rep("black",length(num_list)))+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
#         axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
#         legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
#   ggsave("Converge_Con_0.02_fupdate.pdf",width=10,height=6)
###########################################################################################
#############################Raw Varaiance##################################################
###########################################################################################
load("Nondyn_raw_0.02.RData")

mse_GRNs_raw<-vector()
par_GRNs_raw<-vector()
for(args in 1:93){
  run_id<-rownames(Nondyn_bio_V0.01_B5T5_rep1000_0.02)[as.numeric(args)]
  name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
  for(GRN_ind in 1:8){
    sim_res_GRS<-get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))
    min_ind<-which(unlist(sim_res_GRS[seq(2,1800,6)])==min(unlist(sim_res_GRS[seq(2,1800,6)])))
    mse_GRNs_raw<-rbind(mse_GRNs_raw,unlist(sim_res_GRS[min_ind*6-4]))
    par_GRNs_raw[length(par_GRNs_raw)+1]<-list(unlist(get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))[min_ind*6-5]))
  }
  print(args)
}
sim_GRS_raw<-t(matrix(mse_GRNs_raw,nrow=8))
####Identify select model
sim_GRS_min_raw<-matrix(0,ncol = 8,nrow=93)
for(i in 1:93){
  sim_GRS_min_raw[i,which(sim_GRS_raw[i,]==min(sim_GRS_raw[i,]))]<-1
}
####Estimated GRS
min_list_raw<-apply(sim_GRS_raw,1,function(x) which(x==min(x)))
####ground truth index
min_list_GT<-c(rep(1,7),rep(2,26),rep(3,12),rep(4,8),rep(5,12),rep(6,8),rep(7,12),rep(8,8))
####Likelihood ratio
sim_GRS_F_raw<-rep(0,93)
for(i in 1:93){
  sim_GRS_F_raw[i]<-log2(exp(1))*(min(sim_GRS_raw[i,-min_list_GT[i]])-sim_GRS_raw[i,min_list_GT[i]])
}
#################################################################
################Converge Test####################################
#################################################################
num_list<-seq(100,300,20)
GRS_F_list<-vector()
for(args in 1:93){
  run_id<-rownames(Nondyn_bio_V0.01_B5T5_rep1000_0.02)[as.numeric(args)]
  name_id<-paste(unlist(strsplit(run_id," "))[1],unlist(strsplit(run_id," "))[2],sep="_")
  
  GRS_F_itm<-vector()
  for(i in 1:length(num_list)){
    mse_GRNs<-vector()
    for(GRN_ind in 1:8){
      sim_res_GRS<-get(paste(name_id,"_nondyn_GRN",GRN_ind,sep=""))
      min_ind<-which(unlist(sim_res_GRS[seq(2,6*num_list[i],6)])==min(unlist(sim_res_GRS[seq(2,6*num_list[i],6)])))
      mse_GRNs<-c(mse_GRNs,unlist(sim_res_GRS[min_ind*6-4]))
    }
    GRS_F_itm<-cbind(GRS_F_itm,log2(exp(1))*(min(mse_GRNs[-min_list_GT[args]])-mse_GRNs[min_list_GT[args]]))
  }
  
  GRS_F_list<-rbind(GRS_F_list,GRS_F_itm)
}

Converge_plot<-data.frame(Ite_num=rep(num_list,each=93),Converge=t(matrix(GRS_F_list,nrow=1)))
colnames(Converge_plot)[1:2]<-c("Ite_num","Converge")
Converge_plot$Ite_num<-factor(Converge_plot$Ite_num, levels = num_list)

# pp <- ggplot(Converge_plot, aes(x=Ite_num,y=Converge)) + geom_boxplot(mapping=aes(fill=Ite_num),color=rep("black",length(num_list)),outlier.shape = NA,alpha = 0.3)
# pp+labs(x="",y = "")+
#   scale_fill_manual(values=rep("black",length(num_list)))+
#   geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2))+
#   scale_color_manual(values=rep("black",length(num_list)))+
#   theme_bw()+
#   theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
#         axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
#         legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
#   ggsave("Converge_raw_0.02_fupdate.pdf",width=10,height=6)
##################################################################################################
################################Parametr estimation###############################################
##################################################################################################
load("Nondyn_WT_bio.RData")
Nonredundent_name<-read.table(file = "Nonredundent_name.txt",sep = "\t")
logic_name<-c("1a2a3","1o2o3","1o(2a3)","1a(2o3)","2o(1a3)","2a(1o3)","3o(1a2)","3a(1o2)")
para_name<-paste(expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,1],
                 expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,2],
                 expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,3],sep="")[-27]
GRS_name<-paste(rep(logic_name,each=26),rep(para_name,8))
###########Ground truth
para_list<-as.data.frame(expand.grid(c(0.1,1,10),c(0.1,1,10),c(0.1,1,10)))[-27,]
para_all<-cbind(para_list[rep(c(1:26),8),],1,0.1,0)
para_true<-para_all[which(GRS_name%in%apply(Nonredundent_name,1,as.character)),]
para_true[,4]<-Nondyn_WT_bio
###############TV
para_est_TV<-t(matrix(unlist(par_GRNs_TV[min_list_TV+seq(0,92*8,8)]),nrow=6))
para_fold_TV<-vector()
for(i in 1:93){
  para_fold_TV<-rbind(para_fold_TV,c(10^para_est_TV[i,1:5])/para_true[i,1:5])
} 
para_fold_TV<-cbind(para_fold_TV,para_est_TV[,6])
# save(min_list_TV,para_est_TV,file = "para_est_TV_1228.RData")
###############generate classified prediction
RS_true<-para_true[,1:3]
RS_EST_TV<-10^(para_est_TV[,1:3])
RS_EST_TV[RS_EST_TV<0.5]<-0.1
RS_EST_TV[RS_EST_TV>5]<-10
RS_EST_TV[RS_EST_TV > 0.5 & RS_EST_TV < 5]<-1
###############generate confusion matrix
Confusion_matrix_TV<-matrix(0, ncol=3, nrow=3)
Confusion_matrix_TV[1,1]<-length(which(RS_EST_TV[RS_true==0.1]==0.1))
Confusion_matrix_TV[1,2]<-length(which(RS_EST_TV[RS_true==1]==0.1))
Confusion_matrix_TV[1,3]<-length(which(RS_EST_TV[RS_true==10]==0.1))

Confusion_matrix_TV[2,1]<-length(which(RS_EST_TV[RS_true==0.1]==1))
Confusion_matrix_TV[2,2]<-length(which(RS_EST_TV[RS_true==1]==1))
Confusion_matrix_TV[2,3]<-length(which(RS_EST_TV[RS_true==10]==1))

Confusion_matrix_TV[3,1]<-length(which(RS_EST_TV[RS_true==0.1]==10))
Confusion_matrix_TV[3,2]<-length(which(RS_EST_TV[RS_true==1]==10))
Confusion_matrix_TV[3,3]<-length(which(RS_EST_TV[RS_true==10]==10))

rownames(Confusion_matrix_TV)<-c("Strong","Medium","Weak")
colnames(Confusion_matrix_TV)<-c("Strong","Medium","Weak")
CM_TV<-as.data.frame(as.table(Confusion_matrix_TV))
colnames(CM_TV)[1:2]<-c("Estimated","Ground_Truth")
##############################################################################################
####################################93 GRS##########################################################
##############################################################################################
############################################################
##########################Nonredundant and weak mapping##################################
############################################################
load("Nonredundant_weak_name.RData")
#################
GRS_true<-cbind(min_list_GT,RS_true)
GRS_map<-cbind(GRS_true,c(1:93))
GRS_map<-rbind(GRS_map,c(0,0,0,0,94))
rownames(GRS_map)<-c(as.character(Nonredundent_name$V1),"weak")
#################
add_name_list<-c(weak_name,nonR_1o_name,nonR_1a_name,nonR_2o_name,nonR_2a_name,nonR_3o_name,nonR_3a_name)
logic_add<-match(unlist(lapply(add_name_list,function(x) unlist(strsplit(x," "))[1])),logic_name)
RS_add_raw<-match(unlist(lapply(add_name_list,function(x) unlist(strsplit(x," "))[2])),para_name)
RS_add<-para_list[RS_add_raw,]
GRS_ind_add<-c(rep(94,length(weak_name)),rep(c(32,33),length(nonR_1o_name)/2),c(58,59,64,65,82:85),
               c(28,31,28,31,28,28,28,31,31,31),c(36,39,42,45,82:85),rep(c(16,25),each=5),
               c(36,39,58,59,42,45,64,65))  

GRS_add<-cbind(logic_add,RS_add,GRS_ind_add)
rownames(GRS_add)<-add_name_list
colnames(GRS_map)<-c("logics","Kd1","Kd2","Kd3","GRS_ind")
colnames(GRS_add)<-c("logics","Kd1","Kd2","Kd3","GRS_ind")
GRS_map<-rbind(GRS_map,GRS_add)
##########################################################
##########################################################
##########################################################
load("hclust7_single.RData")
min_list_GT<-c(rep(1,7),rep(2,26),rep(3,12),rep(4,8),rep(5,12),rep(6,8),rep(7,12),rep(8,8))
GRS_est_TV<-cbind(min_list_TV,RS_EST_TV)

Confusion_matrix_NR_GRS_TV<-matrix(0, ncol=93, nrow=94)
for(i in 1:93){
  Confusion_matrix_NR_GRS_TV[GRS_map[row.match(GRS_est_TV[i,],GRS_map[,1:4]),5],i]<-1
}

rownames(Confusion_matrix_NR_GRS_TV)<-c(as.character(Nonredundent_name$V1),"weak")
colnames(Confusion_matrix_NR_GRS_TV)<-Nonredundent_name$V1
CM_NR_GRS_TV<-as.data.frame(as.table(Confusion_matrix_NR_GRS_TV))
colnames(CM_NR_GRS_TV)[1:2]<-c("Estimated","Ground_Truth")
# p <-ggplot(data =CM_NR_GRS_TV,aes(x = Ground_Truth, y = Estimated)) +
#   geom_tile(aes(fill = log(Freq+1)), colour = "white") +
#   labs(x ="Ground Truth", y = "Estimated",fill="Counts")+
#   scale_fill_gradient(low = "white", high = "steelblue") +
#   geom_text(aes(x = Ground_Truth, y = Estimated, label = Freq),size=rel(1.8)) +
#   theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
#         legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
#   ggsave("Confusion_matrix_NR_GRS_0.02_TV_fupdate.pdf",width=10,height=10) 
CM_NR_GRS_TV_ClusterOrder<-as.data.frame(as.table(Confusion_matrix_NR_GRS_TV[c(hclust7_single$order,94),hclust7_single$order]))
colnames(CM_NR_GRS_TV_ClusterOrder)[1:2]<-c("Estimated","Ground_Truth")
# p <-ggplot(data=CM_NR_GRS_TV_ClusterOrder,aes(x = Ground_Truth, y = Estimated)) +
#   geom_tile(aes(fill = log(Freq+1)), colour = "white") +
#   labs(x ="Ground Truth", y = "Estimated",fill="Counts")+
#   scale_fill_gradient(low = "white", high = "steelblue") +
#   geom_text(aes(x = Ground_Truth, y = Estimated, label = Freq),size=rel(1.8)) +
#   theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
#         legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
#   ggsave("Confusion_matrix_NR_GRS_0.02_TV_cluster_order_fupdate.pdf",width=10,height=10) 
############################################################
###############Con##########################################
############################################################
para_est_Con<-t(matrix(unlist(par_GRNs_Con[min_list_Con+seq(0,92*8,8)]),nrow=6))
para_fold_Con<-vector()
for(i in 1:93){
  para_fold_Con<-rbind(para_fold_Con,c(10^para_est_Con[i,1:5])/para_true[i,1:5])
} 
para_fold_Con<-cbind(para_fold_Con,para_est_Con[,6])

RS_EST_Con<-10^(para_est_Con[,1:3])
RS_EST_Con[RS_EST_Con<0.5]<-0.1
RS_EST_Con[RS_EST_Con>5]<-10
RS_EST_Con[RS_EST_Con > 0.5 & RS_EST_Con < 5]<-1

Confusion_matrix_Con<-matrix(0, ncol=3, nrow=3)
Confusion_matrix_Con[1,1]<-length(which(RS_EST_Con[RS_true==0.1]==0.1))
Confusion_matrix_Con[1,2]<-length(which(RS_EST_Con[RS_true==1]==0.1))
Confusion_matrix_Con[1,3]<-length(which(RS_EST_Con[RS_true==10]==0.1))

Confusion_matrix_Con[2,1]<-length(which(RS_EST_Con[RS_true==0.1]==1))
Confusion_matrix_Con[2,2]<-length(which(RS_EST_Con[RS_true==1]==1))
Confusion_matrix_Con[2,3]<-length(which(RS_EST_Con[RS_true==10]==1))

Confusion_matrix_Con[3,1]<-length(which(RS_EST_Con[RS_true==0.1]==10))
Confusion_matrix_Con[3,2]<-length(which(RS_EST_Con[RS_true==1]==10))
Confusion_matrix_Con[3,3]<-length(which(RS_EST_Con[RS_true==10]==10))

rownames(Confusion_matrix_Con)<-c("Strong","Medium","Weak")
colnames(Confusion_matrix_Con)<-c("Strong","Medium","Weak")
CM_Con<-as.data.frame(as.table(Confusion_matrix_Con))
colnames(CM_Con)[1:2]<-c("Estimated","Ground_Truth")
# ggplotConfusionMatrix(CM_Con,"Confusion_matrix_0.02_Con_fupdate","grey50")
#####################################################################
################################GRS Confusion matrix Conventional#####################################
#####################################################################
###################################################################################
###################################################################################
GRS_est_Con<-cbind(min_list_Con,RS_EST_Con)

Confusion_matrix_NR_GRS_Con<-matrix(0, ncol=93, nrow=94)
for(i in 1:93){
  Confusion_matrix_NR_GRS_Con[GRS_map[row.match(GRS_est_Con[i,],GRS_map[,1:4]),5],i]<-1
}

rownames(Confusion_matrix_NR_GRS_Con)<-c(as.character(Nonredundent_name$V1),"weak")
colnames(Confusion_matrix_NR_GRS_Con)<-Nonredundent_name$V1
CM_NR_GRS_Con<-as.data.frame(as.table(Confusion_matrix_NR_GRS_Con))
colnames(CM_NR_GRS_Con)[1:2]<-c("Estimated","Ground_Truth")

# p <-ggplot(data =CM_NR_GRS_Con,aes(x = Ground_Truth, y = Estimated)) +
#   geom_tile(aes(fill = log(Freq+1)), colour = "white") +
#   labs(x ="Ground Truth", y = "Estimated",fill="Counts")+
#   scale_fill_gradient(low = "white", high = "steelblue") +
#   geom_text(aes(x = Ground_Truth, y = Estimated, label = Freq),size=rel(1.8)) +
#   theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
#         legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
#   ggsave("Confusion_matrix_NR_GRS_0.02_Con_fupdate.pdf",width=10,height=10)
##########################################
CM_NR_GRS_Con_ClusterOrder<-as.data.frame(as.table(Confusion_matrix_NR_GRS_Con[c(hclust7_single$order,94),hclust7_single$order]))
colnames(CM_NR_GRS_Con_ClusterOrder)[1:2]<-c("Estimated","Ground_Truth")

# p <-ggplot(data=CM_NR_GRS_Con_ClusterOrder,aes(x = Ground_Truth, y = Estimated)) +
#   geom_tile(aes(fill = log(Freq+1)), colour = "white") +
#   labs(x ="Ground Truth", y = "Estimated",fill="Counts")+
#   scale_fill_gradient(low = "white", high = "steelblue") +
#   geom_text(aes(x = Ground_Truth, y = Estimated, label = Freq),size=rel(1.8)) +
#   theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
#         legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
#   ggsave("Confusion_matrix_NR_GRS_0.02_Con_cluster_order_fupdate.pdf",width=10,height=10)
############################################################
###############raw##########################################
############################################################
para_est_raw<-t(matrix(unlist(par_GRNs_raw[min_list_raw+seq(0,92*8,8)]),nrow=6))
para_fold_raw<-vector()
for(i in 1:93){
  para_fold_raw<-rbind(para_fold_raw,c(10^para_est_raw[i,1:5])/para_true[i,1:5])
} 
para_fold_raw<-cbind(para_fold_raw,para_est_raw[,6])

RS_EST_raw<-10^(para_est_raw[,1:3])
RS_EST_raw[RS_EST_raw<0.5]<-0.1
RS_EST_raw[RS_EST_raw>5]<-10
RS_EST_raw[RS_EST_raw > 0.5 & RS_EST_raw < 5]<-1

Confusion_matrix_raw<-matrix(0, ncol=3, nrow=3)
Confusion_matrix_raw[1,1]<-length(which(RS_EST_raw[RS_true==0.1]==0.1))
Confusion_matrix_raw[1,2]<-length(which(RS_EST_raw[RS_true==1]==0.1))
Confusion_matrix_raw[1,3]<-length(which(RS_EST_raw[RS_true==10]==0.1))

Confusion_matrix_raw[2,1]<-length(which(RS_EST_raw[RS_true==0.1]==1))
Confusion_matrix_raw[2,2]<-length(which(RS_EST_raw[RS_true==1]==1))
Confusion_matrix_raw[2,3]<-length(which(RS_EST_raw[RS_true==10]==1))

Confusion_matrix_raw[3,1]<-length(which(RS_EST_raw[RS_true==0.1]==10))
Confusion_matrix_raw[3,2]<-length(which(RS_EST_raw[RS_true==1]==10))
Confusion_matrix_raw[3,3]<-length(which(RS_EST_raw[RS_true==10]==10))

rownames(Confusion_matrix_raw)<-c("Strong","Medium","Weak")
colnames(Confusion_matrix_raw)<-c("Strong","Medium","Weak")
CM_raw<-as.data.frame(as.table(Confusion_matrix_raw))
colnames(CM_raw)[1:2]<-c("Estimated","Ground_Truth")

ggplotConfusionMatrix(CM_raw,"Confusion_matrix_0.02_raw_fupdate","deepskyblue3")
#####################################################################
################################GRS Confusion matrix Conventional#####################################
#####################################################################
###################################################################################
###################################################################################
load("hclust7_single.RData")
min_list_GT<-c(rep(1,7),rep(2,26),rep(3,12),rep(4,8),rep(5,12),rep(6,8),rep(7,12),rep(8,8))
GRS_est_raw<-cbind(min_list_raw,RS_EST_raw)

Confusion_matrix_NR_GRS_raw<-matrix(0, ncol=93, nrow=94)
for(i in 1:93){
  Confusion_matrix_NR_GRS_raw[GRS_map[row.match(GRS_est_raw[i,],GRS_map[,1:4]),5],i]<-1
}

rownames(Confusion_matrix_NR_GRS_raw)<-c(as.character(Nonredundent_name$V1),"weak")
colnames(Confusion_matrix_NR_GRS_raw)<-Nonredundent_name$V1
CM_NR_GRS_raw<-as.data.frame(as.table(Confusion_matrix_NR_GRS_raw))
colnames(CM_NR_GRS_raw)[1:2]<-c("Estimated","Ground_Truth")

p <-ggplot(data =CM_NR_GRS_raw,aes(x = Ground_Truth, y = Estimated)) +
  geom_tile(aes(fill = log(Freq+1)), colour = "white") +
  labs(x ="Ground Truth", y = "Estimated",fill="Counts")+
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(x = Ground_Truth, y = Estimated, label = Freq),size=rel(1.8)) +
  theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
        legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Confusion_matrix_NR_GRS_0.02_raw_fupdate.pdf",width=10,height=10) 

CM_NR_GRS_raw_ClusterOrder<-as.data.frame(as.table(Confusion_matrix_NR_GRS_raw[c(hclust7_single$order,94),hclust7_single$order]))
colnames(CM_NR_GRS_raw_ClusterOrder)[1:2]<-c("Estimated","Ground_Truth")

p <-ggplot(data=CM_NR_GRS_raw_ClusterOrder,aes(x = Ground_Truth, y = Estimated)) +
  geom_tile(aes(fill = log(Freq+1)), colour = "white") +
  labs(x ="Ground Truth", y = "Estimated",fill="Counts")+
  scale_fill_gradient(low = "white", high = "steelblue") +
  geom_text(aes(x = Ground_Truth, y = Estimated, label = Freq),size=rel(1.8)) +
  theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
        legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Confusion_matrix_NR_GRS_0.02_raw_cluster_order_fupdate.pdf",width=10,height=10) 
##################################################################################
##################################################################################
save.image("Optimization_analysis_0.02_fupdate.RData")
# load("Optimization_analysis_0.02_fupdate.RData")
##################################################################################
#################################################Final Figure#################################
##################################################################################
##################################################################################
##################################################################################
# merge_TV<-Confusion_matrix_NR_GRS_TV[c(hclust7_single$order,94),hclust7_single$order]
# merge_Con<-Confusion_matrix_NR_GRS_Con[c(hclust7_single$order,94),hclust7_single$order]
# merge_raw<-Confusion_matrix_NR_GRS_raw[c(hclust7_single$order,94),hclust7_single$order]

# merge_TV[merge_TV==1]<-2
# marge_both<-merge_TV+merge_Con
# 
# CM_both_ClusterOrder<-as.data.frame(as.table(marge_both))
# colnames(CM_both_ClusterOrder)[1:2]<-c("Estimated","Ground_Truth")
# 
# p <-ggplot(data=CM_both_ClusterOrder,aes(x = Ground_Truth, y = Estimated)) +
#   geom_tile(aes(fill = factor(Freq))) +
#   labs(x ="Ground Truth", y = "Estimated")+
#   scale_fill_manual(values=c("white","black","red","#8B0000")) +
#   theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
#         legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
#   ggsave("Confusion_matrix_GRS_0.02_TV_cluster_order_merge_fupdate.pdf",width=10,height=10) 
##################################################################################
##################################################################################
#################################################Final Figure#################################
##################################################################################
##################################################################################
merge_TV<-Confusion_matrix_NR_GRS_TV[c(hclust7_single$order,94),hclust7_single$order]
merge_Con<-Confusion_matrix_NR_GRS_Con[c(hclust7_single$order,94),hclust7_single$order]
merge_raw<-Confusion_matrix_NR_GRS_raw[c(hclust7_single$order,94),hclust7_single$order]

merge_Con[merge_Con==1]<-3
merge_raw[merge_raw==1]<-5
marge_all<-merge_TV+merge_Con+merge_raw

CM_both_ClusterOrder<-as.data.frame(as.table(marge_all))
colnames(CM_both_ClusterOrder)[1:2]<-c("Estimated","Ground_Truth")

p <-ggplot(data=CM_both_ClusterOrder,aes(x = Ground_Truth, y = Estimated)) +
  geom_tile(aes(fill = factor(Freq))) +
  labs(x ="Ground Truth", y = "Estimated")+
  scale_fill_manual(values=c("white","red","black","red4","deepskyblue3","navy")) +
  theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
        legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Confusion_matrix_GRS_0.02_mergeAll_fupdate.pdf",width=10,height=10)
##################################################################################
merge_org_TV<-Confusion_matrix_NR_GRS_TV
merge_org_Con<-Confusion_matrix_NR_GRS_Con
merge_org_raw<-Confusion_matrix_NR_GRS_raw

merge_org_Con[merge_org_Con==1]<-3
merge_org_raw[merge_org_raw==1]<-5
merge_org_all<-merge_org_TV+merge_org_Con+merge_org_raw

CM_both_ClusterOrder<-as.data.frame(as.table(merge_org_all))
colnames(CM_both_ClusterOrder)[1:2]<-c("Estimated","Ground_Truth")

p <-ggplot(data=CM_both_ClusterOrder,aes(x = Ground_Truth, y = Estimated)) +
  geom_tile(aes(fill = factor(Freq))) +
  labs(x ="Ground Truth", y = "Estimated")+
  scale_fill_manual(values=c("white","red","black","red4","deepskyblue3","navy")) +
  theme(axis.text.x=element_text(angle = 90,vjust = 0.5,face="bold",size=rel(0.5)),axis.text.y=element_text(size=rel(0.5),face="bold"),axis.title=element_text(size=rel(1),face="bold"),
        legend.position="none",plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Confusion_matrix_GRS_0.02_mergeAll_original_fupdate.pdf",width=10,height=10)
##################################################################################
##################################################################################
##################################################################################
para_est_plot<-data.frame(para=rep(c("Kd1_TV","Kd2_TV","Kd3_TV","ksyn_TV","kdeg_TV","Kd1_Con","Kd2_Con","Kd3_Con","ksyn_Con","kdeg_Con","Kd1_raw","Kd2_raw","Kd3_raw","ksyn_raw","kdeg_raw"),each=93),est_rel=c(para_fold_TV[,1],para_fold_TV[,2],para_fold_TV[,3],para_fold_TV[,4],para_fold_TV[,5],para_fold_Con[,1],para_fold_Con[,2],para_fold_Con[,3],para_fold_Con[,4],para_fold_Con[,5],para_fold_raw[,1],para_fold_raw[,2],para_fold_raw[,3],para_fold_raw[,4],para_fold_raw[,5]))
para_est_plot$para<-factor(para_est_plot$para, levels = c("Kd1_TV","Kd1_Con","Kd1_raw","Kd2_TV","Kd2_Con","Kd2_raw","Kd3_TV","Kd3_Con","Kd3_raw","ksyn_TV","ksyn_Con","ksyn_raw","kdeg_TV","kdeg_Con","kdeg_raw"))

per_cal<-function(x){
  length(intersect(which(x<2),which(x>0.5)))/93
}

c(per_cal(para_fold_TV[,1]),per_cal(para_fold_TV[,2]),per_cal(para_fold_TV[,3]),per_cal(para_fold_Con[,1]),per_cal(para_fold_Con[,2]),per_cal(para_fold_Con[,3]),per_cal(para_fold_raw[,1]),per_cal(para_fold_raw[,2]),per_cal(para_fold_raw[,3]))

# para_est_plot<-data.frame(para=rep(c("Kd1","Kd2","Kd3","ksyn","kdeg","Kd1","Kd2","Kd3","ksyn","kdeg","Kd1","Kd2","Kd3","ksyn","kdeg"),each=93),
#                           model = rep(c("TV","Con","raw"), each=93*5),est_rel=c(para_fold_TV[,1],para_fold_TV[,2],para_fold_TV[,3],para_fold_TV[,4],para_fold_TV[,5],para_fold_Con[,1],para_fold_Con[,2],para_fold_Con[,3],para_fold_Con[,4],para_fold_Con[,5],para_fold_raw[,1],para_fold_raw[,2],para_fold_raw[,3],para_fold_raw[,4],para_fold_raw[,5]))
# para_est_plot$model<-factor(para_est_plot$model, levels = c("TV","Con","raw"))

#########################################Plot Para estimation#########################################################
# pp <- ggplot(para_est_plot, aes(x=para,y=est_rel, fill=model,col=model)) + geom_boxplot(outlier.shape = NA,alpha = 0.3)
# pp <- pp +  geom_hline(yintercept=1,color = adjustcolor("orange", alpha.f = 0.5), size=1,linetype = "dashed") 
# pp <- pp+ scale_y_continuous(trans='log10', limits=c(0.001, 1000)) + scale_fill_manual(values=rep(c("red","black","deepskyblue3"),5))+
#   geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2))+
#   scale_color_manual(values=c("red","black","deepskyblue3"))
#   #scale_x_discrete(labels=rep(c("Kd1","Kd2","Kd3","ksyn","kdeg")))
# pp <- pp+theme_bw()+
#   theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
#         axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
#         legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))
# pp

pp <- ggplot(para_est_plot, aes(x=para,y=est_rel)) + geom_boxplot(mapping=aes(fill=para),outlier.shape = NA,color=rep(c("red","black","deepskyblue3"),5),alpha = 0.3)
pp+labs(x="",y = "")+
  geom_hline(yintercept=1,color = adjustcolor("orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
  scale_y_continuous(trans='log10', limits=c(0.001, 1000))+
  scale_fill_manual(values=rep(c("red","black","deepskyblue3"),5))+
  geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),aes(color = para))+
  scale_color_manual(values=rep(c("red","black","deepskyblue3"),5))+
  scale_x_discrete(labels=rep(c("Kd1","Kd2","Kd3","ksyn","kdeg"),each=3))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Para_est_0.02_fupdate.pdf",width=8,height=6)
##################################################################################
#####################################Kds#############################################
##################################################################################
para_est_plot_Kds<-data.frame(para=rep(c("Kd1_TV","Kd2_TV","Kd3_TV","Kd1_Con","Kd2_Con","Kd3_Con","Kd1_raw","Kd2_raw","Kd3_raw"),each=93),est_rel=c(para_fold_TV[,1],para_fold_TV[,2],para_fold_TV[,3],para_fold_Con[,1],para_fold_Con[,2],para_fold_Con[,3],para_fold_raw[,1],para_fold_raw[,2],para_fold_raw[,3]))
para_est_plot_Kds$para<-factor(para_est_plot_Kds$para, levels = c("Kd1_TV","Kd1_Con","Kd1_raw","Kd2_TV","Kd2_Con","Kd2_raw","Kd3_TV","Kd3_Con","Kd3_raw"))
#########################################
pp <- ggplot(para_est_plot_Kds, aes(x=para,y=est_rel)) + geom_boxplot(mapping=aes(fill=para),outlier.shape = NA,color=rep(c("red","black","deepskyblue3"),3),alpha = 0.3)
pp+labs(x="",y = "")+
  geom_hline(yintercept=1,color = adjustcolor("orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
  scale_y_continuous(trans='log10')+
  scale_fill_manual(values=rep(c("red","black","deepskyblue3"),3))+
  geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),aes(color = para))+
  scale_color_manual(values=rep(c("red","black","deepskyblue3"),3))+
  scale_x_discrete(labels=rep(c("Kd1","Kd2","Kd3"),each=3))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Para_est_Kds_0.02_fupdate.pdf",width=6,height=6)
############################basal parameters
para_est_k0_plot<-data.frame(para=rep(c("k0_TV","k0_Con","k0_raw"),each=93),est_rel=c(para_fold_TV[,6],para_fold_Con[,6],para_fold_raw[,6]))
para_est_k0_plot$para<-factor(para_est_k0_plot$para, levels = c("k0_TV","k0_Con","k0_raw"))

pp <- ggplot(para_est_k0_plot, aes(x=para,y=est_rel)) + geom_boxplot(mapping=aes(fill=para),outlier.shape = NA,color=c("red","black","deepskyblue3"),alpha = 0.3)
pp+labs(x="",y = "")+
  geom_hline(yintercept=0,color = adjustcolor("orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
  scale_fill_manual(values=c("red","black","deepskyblue3"))+
  geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),aes(color = para))+
  scale_color_manual(values=c("red","black","deepskyblue3"))+
  scale_x_discrete(labels=rep(c("k0","k0"),each=3))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Para_est_k0_0.02_fupdate.pdf",width=4,height=6)

pp <- ggplot(para_est_k0_plot, aes(x=para,y=est_rel)) + geom_boxplot(mapping=aes(fill=para),outlier.shape = NA,color=c("red","black","deepskyblue3"),alpha = 0.3)
pp+labs(x="",y = "")+
  geom_hline(yintercept=0,color = adjustcolor("orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
  scale_fill_manual(values=c("red","black","deepskyblue3"))+
  geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),aes(color = para))+
  scale_color_manual(values=c("red","black","deepskyblue3"))+
  scale_x_discrete(labels=rep(c("k0","k0"),each=3))+
  ylim(0,0.035)+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Para_est_k0_0.02_break_fupdate.pdf",width=4,height=6)
########################################################################################
########################################################################################
########################################################################################
Iden_plot<-data.frame(error_model=rep(c("Con","TV","raw"),each=93),Iden=c(sim_GRS_F_Con,sim_GRS_F_TV,sim_GRS_F_raw))

pp <- ggplot(Iden_plot, aes(x=error_model,y=Iden)) + geom_boxplot(mapping=aes(fill=error_model),outlier.shape = NA,color=c("black","red","deepskyblue3"),alpha = 0.3)
pp+labs(x="",y = "")+
  geom_hline(yintercept=0,color = adjustcolor( "orange", alpha.f = 0.5), size=1,linetype = "dashed") + 
  scale_fill_manual(values=c("black","red","deepskyblue3"))+
  geom_jitter(shape=16, position=position_jitter(height = 0, width = 0.2),aes(color = error_model))+
  scale_color_manual(values=c("black","red","deepskyblue3"))+
  theme_bw()+
  theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))  +
  ggsave("Identifiability_0.02_fupdate.pdf",width=4,height=6)
########################################################################################
########################################Bar plot################################################
########################################################################################
########################################################################################
Iden_number<-data.frame(error_model=c("TV","Con","raw"),Iden_n=c(93,93,92))
Iden_number$error_model<-factor(Iden_number$error_model, levels = c("TV","Con","raw"))

pp <- ggplot(Iden_number, aes(x=error_model,y=Iden_n, fill=error_model,col=error_model)) +geom_bar(stat="identity")
pp <- pp + geom_hline(yintercept=93,color = adjustcolor("orange", alpha.f = 0.5), size=1,linetype = "dashed")
pp <- pp+ scale_y_continuous(breaks=c(50, 93),limits=c(0,100)) + scale_color_manual(values=c("red","black","deepskyblue3"), guide = F)+ scale_fill_manual(values=c("red","black","deepskyblue3"))+labs(x="",y = "")
pp <- pp+theme_bw()+
  theme(axis.text.x=element_text(angle = 90, hjust = 1,vjust=0.5,face="bold",size=rel(1.5)),axis.text.y=element_text(size=rel(1.5),face="bold"),
        axis.title=element_text(size=rel(2),face="bold"),legend.text=element_text(size=rel(1.5)),
        legend.title =element_text(size=rel(1.5)),plot.margin = unit(c(1,1,1,1), "cm"))
pp+ggsave("Iden_gate_0.02_fupdate.pdf",width=4.5,height=6)
########################################################################################
########################################################################################
##########################################Plot the simulation###########################
########################################################################################
# rm(list=ls())
library("MASS")
library(deSolve)
library(signal)
library(pheatmap)
source("Nondyn_opt_newbasalnew.R")
load("Nondyn_bio_V0.01_B5T5_rep1000_0.02.RData")
load("para_est_TV_1228.RData")
load("Nondyn_active_GT.RData")
# load("Prior_TV_bio_V0.01_B5T5_rep1000_0.02_fupdate.RData")
# var_est<-Prior_TV_bio_V0.01_B5T5_rep1000_0.02_fupdate$var_est_n
# var_est_g<-Prior_TV_bio_V0.01_B5T5_rep1000_0.02_fupdate$var_est_g
Nondyn_V0.01_TB5_TT5_rep1000_p<-cbind(Nondyn_bio_V0.01_B5T5_rep1000_0.02[,1:104],Nondyn_bio_V0.01_B5T5_rep1000_0.02[,105:208])
Nondyn_V0.01_TB5_TT5_Average<-(Nondyn_V0.01_TB5_TT5_rep1000_p[,1:104]+Nondyn_V0.01_TB5_TT5_rep1000_p[,105:208])/2
#####################################################
#####################################################True TF activities
Basal_level=0
Medium_strength<-0.5
input_list_raw<-as.data.frame(expand.grid(c(1,Medium_strength,Basal_level),c(1,Medium_strength,Basal_level),c(1,Medium_strength,Basal_level)))[-27,]
####change deltaN
deltaN<-1
set.seed(1)
input_error<-matrix(rnorm(3*26,mean=0,sd=0.707*deltaN*0.02),nrow=26) 
input_list<-input_list_raw*(1+input_error)
#####################################################
#####################################################
#####################################################Function test
#####################################################
#############Error for RNA
####likelihood of RNA dynamics fitting
caRNA_time<-c(0,15,30,60)
############################################################
dynamic_exp_sim<-vector()

for(para_ind in 1:93){
  Kd1<-10^(para_est_TV[para_ind,1])
  Kd2<-10^(para_est_TV[para_ind,2])
  Kd3<-10^(para_est_TV[para_ind,3])
  ksyn<-10^(para_est_TV[para_ind,4])
  kdeg<-10^(para_est_TV[para_ind,5])
  k0<-para_est_TV[para_ind,6]
  GRN_ind<-min_list_TV[para_ind]
  
  parS<-c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,k0=k0,ksyn=ksyn,kdeg=kdeg)
  R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg
  Gene_exp_sim<-c(get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[1,1],h2=input_list[1,2],h3=input_list[1,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[2,1],h2=input_list[2,2],h3=input_list[2,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[3,1],h2=input_list[3,2],h3=input_list[3,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[4,1],h2=input_list[4,2],h3=input_list[4,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[5,1],h2=input_list[5,2],h3=input_list[5,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[6,1],h2=input_list[6,2],h3=input_list[6,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[7,1],h2=input_list[7,2],h3=input_list[7,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[8,1],h2=input_list[8,2],h3=input_list[8,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[9,1],h2=input_list[9,2],h3=input_list[9,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[10,1],h2=input_list[10,2],h3=input_list[10,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[11,1],h2=input_list[11,2],h3=input_list[11,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[12,1],h2=input_list[12,2],h3=input_list[12,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[13,1],h2=input_list[13,2],h3=input_list[13,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[14,1],h2=input_list[14,2],h3=input_list[14,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[15,1],h2=input_list[15,2],h3=input_list[15,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[16,1],h2=input_list[16,2],h3=input_list[16,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[17,1],h2=input_list[17,2],h3=input_list[17,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[18,1],h2=input_list[18,2],h3=input_list[18,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[19,1],h2=input_list[19,2],h3=input_list[19,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[20,1],h2=input_list[20,2],h3=input_list[20,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[21,1],h2=input_list[21,2],h3=input_list[21,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[22,1],h2=input_list[22,2],h3=input_list[22,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[23,1],h2=input_list[23,2],h3=input_list[23,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[24,1],h2=input_list[24,2],h3=input_list[24,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[25,1],h2=input_list[25,2],h3=input_list[25,3])),
                  get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[26,1],h2=input_list[26,2],h3=input_list[26,3])))
  dynamic_exp_sim<-rbind(dynamic_exp_sim,Gene_exp_sim)
}
############################################################
############################################################
pdf("Nondyn_fitting_performance_0.02_1008.pdf",width=10,height=6)
pheatmap(cbind(Nondyn_V0.01_TB5_TT5_rep1000_p,dynamic_exp_sim),gaps_col = c(104,208),cellheight=4,cellwidth=2,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()
############################################################
############################################################
pdf("Nondyn_fitting_performance_0.02_full_1228.pdf",width=7,height=6)
pheatmap(cbind(Nondyn_active_GT,Nondyn_V0.01_TB5_TT5_rep1000_p,dynamic_exp_sim),gaps_col = c(104,208,312),cellheight=4,cellwidth=1,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()