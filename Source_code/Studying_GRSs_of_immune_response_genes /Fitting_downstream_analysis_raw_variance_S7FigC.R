rm(list=ls())
# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Exp_fit")
library(ggplot2)
library(pheatmap)

mse_GRNs<-vector()
par_GRNs<-vector()
for(Gene_ind in c(7,16)){
  if(file.exists(paste("dyn_Gene",Gene_ind,".RData",sep=""))){
    load(paste("dyn_Gene",Gene_ind,".RData",sep=""))
    for(GRN_ind in 1:8){
      sim_res_GRS<-get(paste("dyn_Gene",Gene_ind,"_Logic",GRN_ind,sep=""))

      min_ind<-which(unlist(sim_res_GRS[seq(2,150,5)])==min(unlist(sim_res_GRS[seq(2,150,5)])))[1]
      mse_GRNs<-rbind(mse_GRNs,unlist(sim_res_GRS[min_ind*5-3]))
      par_GRNs[length(par_GRNs)+1]<-list(unlist(sim_res_GRS[min_ind*5-4]))
    }
    print(Gene_ind)
  }
}

GRS_fit<-t(matrix(mse_GRNs,nrow=8))
GRS_best_fit<-apply(GRS_fit,1,min)
############################################################
############################################################
####################Generate fitted gene expression data####
############################################################
############################################################
library(deSolve)
source("Dyn_fit_equa_V3.R")
load("Experimental_data.RData")
int_length<-10^-3
Basal_level<-0
#############Error for RNA
load("WT_LipidA_GeneExp.Prior.RData")
load("MAPKi_LipidA_GeneExp.Prior.RData")
load("MYD88_LipidA_GeneExp.Prior.RData")
#############Using WT condition to estimate TRIF KO and Pam3csk4, whihc doesn't have 2 replicates  
var_est_g_list<-rbind(WT_LipidA_GeneExp.Prior$var_est_g,WT_LipidA_GeneExp.Prior$var_est_g,MYD88_LipidA_GeneExp.Prior$var_est_g,MAPKi_LipidA_GeneExp.Prior$var_est_g,WT_LipidA_GeneExp.Prior$var_est_g)
var_est_n_list<-rbind(WT_LipidA_GeneExp.Prior$var_est_n,WT_LipidA_GeneExp.Prior$var_est_n,MYD88_LipidA_GeneExp.Prior$var_est_n,MAPKi_LipidA_GeneExp.Prior$var_est_n,WT_LipidA_GeneExp.Prior$var_est_n)
###############################calculate slope#################
slope_ind<-round(2*(c(rep(0,5),var_est_n_list[,4]))^0.5)
slope_ind[slope_ind>60]<-60
slope_ind[slope_ind<1]<-1
num_cond<-dim(Exp_data_mean)[2]/4

caRNA_time_org<-c(c(0,15,30,60),as.vector(t(rbind(c(15,15,30,30,60,60),c(15,15,30,30,60,60))[rep(1,num_cond*2),]+c(-slope_ind,slope_ind))))
caRNA_time_org[caRNA_time_org<0]<-0
caRNA_time_org[caRNA_time_org>80]<-80
slope_interval_g<-(caRNA_time_org[seq(6,length(caRNA_time_org),2)]-caRNA_time_org[seq(5,length(caRNA_time_org),2)])[1:(num_cond*3)]
slope_interval_n<-(caRNA_time_org[seq(6,length(caRNA_time_org),2)]-caRNA_time_org[seq(5,length(caRNA_time_org),2)])[(num_cond*3+1):(num_cond*6)]
caRNA_time<-sort(caRNA_time_org)
back_order<-rep(order(order(caRNA_time_org)),num_cond)+rep(seq(0,64*(num_cond-1),64),each=64)
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
####################Generate fitted gene expression data###############
data_fit_list<-vector()
Gene_list_ind<-c(7,16)
for(Gene_ind in 1:2){
  exp_data<-as.numeric(Exp_data_mean[Gene_list_ind[Gene_ind],])
  ##add sudo counts for zero counts
  exp_data[exp_data==0]<-0.1
  
  exp_data_basal<-rep(exp_data[seq(1,20,4)],each=64)
  num_cond<-length(exp_data)/4
  data_fit_list<-rbind(data_fit_list,exp_data)
  
  for(GRN_ind in 1:8){
    Kd1<-as.numeric(10^(par_GRNs[[GRN_ind+(Gene_ind-1)*8]][1]))
    Kd2<-as.numeric(10^(par_GRNs[[GRN_ind+(Gene_ind-1)*8]][2]))
    Kd3<-as.numeric(10^(par_GRNs[[GRN_ind+(Gene_ind-1)*8]][3]))
    ksyn<-as.numeric(10^(par_GRNs[[GRN_ind+(Gene_ind-1)*8]][4]))
    kdeg<-as.numeric(10^(par_GRNs[[GRN_ind+(Gene_ind-1)*8]][5]))
    
    R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,Basal_level))/kdeg
    
    model_data<-matrix(ode(y=c(RNA1=R_basal,RNA2=R_basal,RNA3=R_basal,RNA4=R_basal,RNA5=R_basal), 
                           times = caRNA_time,func =get(paste("P1_GRN",GRN_ind,"_fit",sep="")),parms = c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,ksyn=ksyn,kdeg=kdeg),hmax=1)[,-1],ncol=1)
    
    data_fit<-(model_data+exp_data_basal)[back_order][data_ind]
    data_fit_list<-rbind(data_fit_list,data_fit)
  }
  
}

target_fun_list<-mse_GRNs
############################################################
############################################################
#########Map GRS to 17 Logics################################
############################################################
############################################################
#############Discrete regulation strength into 4 categories###########
RS_matrix<-rbind(c(NA,NA,NA),c(NA,NA,NA))[rep(1,2*8),]

RS_raw<-t(matrix(unlist(par_GRNs),nrow=5))[,1:3]
RS_matrix[which(10^(RS_raw)<=0.5)]<-"S"
RS_matrix[intersect(which(10^(RS_raw)<=5),which(10^(RS_raw)>0.5))]<-"M"
RS_matrix[intersect(which(10^(RS_raw)<=100),which(10^(RS_raw)>5))]<-"W"
RS_matrix[which(10^(RS_raw)>100)]<-"N"
############################################################
############################################################
#########Plot negative log likelihood for all the genes#####
############################################################
############################################################
########################################################################################################################
########################################visualize for two individual genes##############################################
########################################################################################################################
fit_threshold<-100

RS_raw_plot<-RS_raw
RS_raw_plot[RS_raw_plot>=2]<-2
RS_raw_plot[RS_raw_plot<=(-1)]<-(-1)

colfunc <- colorRampPalette(c("white", "firebrick"))
logic_disp<-rbind(c("A","A","A"),c("O","O","O"),c("O","O","A"),c("A","A","O"),c("O","A","O"),c("A","O","A"),c("A","O","O"),c("O","A","A"))[rep(c(1:8),2),]

annotation_row_fit_plot<-as.data.frame(cbind(logic_disp,RS_raw_plot)[,c(6,5,4,3,2,1)])
annotation_row_fit_plot<-cbind(target_fun_list,annotation_row_fit_plot)
annotation_row_fit_plot[,2:4]<-RS_raw_plot[,c(3,2,1)]

annotation_row_plot<-as.data.frame(rbind(c(NA,NA,NA,NA,NA,NA,NA),c(NA,NA,NA,NA,NA,NA,NA))[rep(1,9*2),])
annotation_row_plot[c(1:dim(annotation_row_plot)[1])[-seq(1,dim(annotation_row_plot)[1],9)],]<-annotation_row_fit_plot

colnames(annotation_row_plot)<-c("LL","TF3","TF2","TF1","TF23","TF13","TF12")
rownames(annotation_row_plot)<-paste(rep(rownames(Exp_data_mean)[c(7,16)],each=9),rep(c("",c("  1a2a3","  1o2o3","  1o(2a3)","  1a(2o3)","  2o(1a3)","  2a(1o3)","  3o(1a2)","  3a(1o2)")),2),sep="")

anno_colors <- list("TF3"=c("firebrick","white"),
                    "TF2"=c("firebrick","white"),
                    "TF1"=c("firebrick","white"),
                    "TF23"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF13"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF12"=c("A"="#E69F00","O"="#56B4E9"),
                    "LL"=c("gold2", "grey"))

data_plot<-data_fit_list
data_plot<-rbind(data_plot,"min"=rep(0.122304,20),"max"=rep(73.26923,20))

anno_plot<-annotation_row_plot
anno_plot[which(anno_plot[,1]<fit_threshold),1]<-0
anno_plot[which(anno_plot[,1]>fit_threshold),1]<-1

range_keep<-as.data.frame(rbind("min"=c(annotation_row_plot[2,1],c(-1,-1,-1),c("A","A","A")),"max"=c(annotation_row_plot[2,1],c(2,2,2),c("A","A","A"))))
range_keep[,1:4]<-rbind(c(0,c(-1,-1,-1)),c(0,c(2,2,2)))
colnames(range_keep)<-colnames(anno_plot)
anno_plot<-rbind(anno_plot,range_keep)

rownames(data_plot)[c(2:9,11:18)]<-rownames(anno_plot)[c(2:9,11:18)]
Gene_order<-c(1,order(target_fun_list[1:8])+1,10,order(target_fun_list[9:16])+10,19,20)

pdf("Exp_fitting_two_genes_Srgn_Ppp1r15a_revise_raw_variance.pdf",width=6,height=13)
pheatmap(log2(data_plot[Gene_order,]+0.01),cellheight=8, cellwidth = 3,gaps_row = c(9,18),gaps_col = seq(4,16,4),
         legend=T,
         show_rownames = T,show_colnames = F,cluster_cols = F,cluster_rows = F,font_size=0.1,
         annotation_row = anno_plot,annotation_legend = T, annotation_names_row = F,annotation_colors = anno_colors)
graphics.off()