rm(list=ls())
# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Exp_fit")
library(ggplot2)
library(pheatmap)

# mse_GRNs<-vector()
# par_GRNs<-vector()
# for(Gene_ind in 1:159){
#   if(file.exists(paste("dyn_Gene",Gene_ind,".RData",sep=""))){
#     load(paste("dyn_Gene",Gene_ind,".RData",sep=""))
#     for(GRN_ind in 1:8){
#       sim_res_GRS<-get(paste("dyn_Gene",Gene_ind,"_Logic",GRN_ind,sep=""))
#       
#       min_ind<-which(unlist(sim_res_GRS[seq(2,150,5)])==min(unlist(sim_res_GRS[seq(2,150,5)])))[1]
#       mse_GRNs<-rbind(mse_GRNs,unlist(sim_res_GRS[min_ind*5-3]))
#       par_GRNs[length(par_GRNs)+1]<-list(unlist(sim_res_GRS[min_ind*5-4]))
#     }
#     print(Gene_ind)
#   }
# }
# 
# GRS_fit<-t(matrix(mse_GRNs,nrow=8))
# GRS_best_fit<-apply(GRS_fit,1,min)
# save.image("Data_output_final.RData")
load("Data_output_final.RData")
############################################################
############################################################
####################Generate fitted gene expression data####
############################################################
############################################################
library(deSolve)

source("Dyn_fit_equa_V3.R")

load("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Exp_fit/Experimental_data.RData")
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
######################target function########
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
####################Generate fitted gene expression data###############
target_fun_list<-vector()
data_fit_list<-vector()

for(Gene_ind in 1:159){
  exp_data<-as.numeric(Exp_data_mean[Gene_ind,])
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
    
    Total_mse<-target_fun(matrix(ode(y=c(RNA1=R_basal,RNA2=R_basal,RNA3=R_basal,RNA4=R_basal,RNA5=R_basal), 
                                     times = caRNA_time,func = get(paste("P1_GRN",GRN_ind,"_fit",sep="")),parms = c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,ksyn=ksyn,kdeg=kdeg),hmax=1)[,-1],ncol=1))
    
    target_fun_list<-c(target_fun_list,Total_mse)
    
    data_fit<-(model_data+exp_data_basal)[back_order][data_ind]
    data_fit_list<-rbind(data_fit_list,data_fit)
  }
  
}

target_fun_list_thre<-target_fun_list
############################################################
############################################################
#########Map GRS to 17 Logics################################
############################################################
############################################################
#############Discrete regulation strength into 4 categories###########
RS_matrix<-rbind(c(NA,NA,NA),c(NA,NA,NA))[rep(1,159*8),]

RS_raw<-t(matrix(unlist(par_GRNs),nrow=5))[,1:3]
RS_matrix[which(10^(RS_raw)<=0.5)]<-"S"
RS_matrix[intersect(which(10^(RS_raw)<=5),which(10^(RS_raw)>0.5))]<-"M"
RS_matrix[intersect(which(10^(RS_raw)<=100),which(10^(RS_raw)>5))]<-"W"
RS_matrix[which(10^(RS_raw)>100)]<-"N"
########################################Map back to the 17 GRS###########
Mapped_logics<-rep(c(10:17),159)
for(i in 1:length(Mapped_logics)){
  if(any(RS_matrix[i,]=="N")){
    
    if(RS_matrix[i,2]=="N" & RS_matrix[i,3]=="N" & i%%8==2) { Mapped_logics[i]<- 1 }
    if((RS_matrix[i,2]=="N" | RS_matrix[i,3]=="N") & i%%8==3) { Mapped_logics[i]<- 1 }
    
    if(RS_matrix[i,1]=="N" & RS_matrix[i,3]=="N" & i%%8==2) { Mapped_logics[i]<- 2 }
    if((RS_matrix[i,1]=="N" | RS_matrix[i,3]=="N") & i%%8==5) { Mapped_logics[i]<- 2 }
    
    if(RS_matrix[i,1]=="N" & RS_matrix[i,2]=="N" & i%%8==2) { Mapped_logics[i]<- 3 }
    if((RS_matrix[i,1]=="N" | RS_matrix[i,2]=="N") & i%%8==7) { Mapped_logics[i]<- 3 }
    
    if(RS_matrix[i,3]=="N" & any(i%%8==c(4,6,7))) { Mapped_logics[i]<- 4 }
    if(RS_matrix[i,3]=="N" & i%%8==2) { Mapped_logics[i]<- 5 }
    
    if(RS_matrix[i,1]=="N" & any(i%%8==c(3,6,8))) { Mapped_logics[i]<- 6 }
    if(RS_matrix[i,1]=="N" & i%%8==2) { Mapped_logics[i]<- 7 }
    
    if(RS_matrix[i,2]=="N" & any(i%%8==c(4,5,8))) { Mapped_logics[i]<- 8 }
    if(RS_matrix[i,2]=="N" & i%%8==2) { Mapped_logics[i]<- 9 }
  }
}

Mapped_logics_M<-t(matrix(Mapped_logics,nrow=8))
################################Assigned logics################################################
LL_list<-t(matrix(target_fun_list,nrow=8))
LL_logic_gates<-matrix(NA, ncol=17,nrow=159)
for(i in 1:159){
  LL_logic_gates[i,Mapped_logics_M[i,order(LL_list[i,],decreasing = T)]]<-sort(LL_list[i,],decreasing = T)
}
rownames(LL_logic_gates)<-rownames(Exp_data_mean)

load("hire_order.RData")
cols<- colorRampPalette(c("Yellow", "Black"))(30)
pdf("Exp_fitting_logics_supplementary.pdf",width=5,height=10)
pheatmap(LL_logic_gates[hire_order$tree_row$order,],cellheight=4, cellwidth = 4,gaps_col = c(3,9),
         legend=T,color = cols,na_col = "white",
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()

save.image("fitting_analysis.RData")
# load("fitting_analysis.RData")
############################################################
############################################################
#########Plot negative log likelihood for all the genes#####
############################################################
############################################################
fit_threshold<-100
############################################################
GRS_LL<-as.vector(t(LL_logic_gates[order(apply(LL_logic_gates,1,function(x) min(x,na.rm = T))),]))

pdf("GRS_Fit_LL.pdf",width=5,height=5)
plot(rep(c(1:159),each=17),GRS_LL,xlab="",ylab="",pch=16,col="grey",cex=0.5)
abline(h=fit_threshold,lty=2,col="orange",lwd=2.5)
graphics.off()
############################################################
############################################################
#########Examine of fitted GRS for all the genes#####
############################################################
############################################################
colfunc <- colorRampPalette(c("white", "firebrick"))
logic_disp<-rbind(c("A","A","A"),c("O","O","O"),c("O","O","A"),c("A","A","O"),c("O","A","O"),c("A","O","A"),c("A","O","O"),c("O","A","A"))[rep(c(1:8),159),]

annotation_row_fit<-as.data.frame(cbind(logic_disp,RS_matrix)[,c(6,5,4,3,2,1)])

annotation_row_fit<-cbind(target_fun_list_thre,annotation_row_fit)

annotation_row<-as.data.frame(rbind(c(NA,NA,NA,NA,NA,NA,NA),c(NA,NA,NA,NA,NA,NA,NA))[rep(1,9*159),])
annotation_row[c(1:dim(annotation_row)[1])[-seq(1,dim(annotation_row)[1],9)],]<-annotation_row_fit

colnames(annotation_row)<-c("LL","TF3","TF2","TF1","TF23","TF13","TF12")
rownames(annotation_row)<-paste(rep(rownames(Exp_data_mean),each=9),rep(c("",c("  1a2a3","  1o2o3","  1o(2a3)","  1a(2o3)","  2o(1a3)","  2a(1o3)","  3o(1a2)","  3a(1o2)")),159),sep="")

anno_colors <- list("TF3"=c("S"=colfunc(4)[4],"M"=colfunc(4)[3],"W"=colfunc(4)[2],"N"=colfunc(4)[1]),
                    "TF2"=c("S"=colfunc(4)[4],"M"=colfunc(4)[3],"W"=colfunc(4)[2],"N"=colfunc(4)[1]),
                    "TF1"=c("S"=colfunc(4)[4],"M"=colfunc(4)[3],"W"=colfunc(4)[2],"N"=colfunc(4)[1]),
                    "TF23"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF13"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF12"=c("A"="#E69F00","O"="#56B4E9"),
                    "LL"=c("Yellow", "Black"))

rownames(data_fit_list)<-paste(rep(rownames(Exp_data_mean),each=9),rep(c("",c("  1a2a3","  1o2o3","  1o(2a3)","  1a(2o3)","  2o(1a3)","  2a(1o3)","  3o(1a2)","  3a(1o2)")),159),sep="")
####################reorder of the results
####oreder of inter GRS
GRN_inter_order<-c(1:(159*9))
GRN_inter_order[c(1:(159*9))[-seq(1,159*9,9)]]<-as.vector((apply(LL_list,1,order))+1)+rep(seq(0,158*9,9),each=8)
data_fit_order1<-data_fit_list[GRN_inter_order,]
annotation_row_order1<-annotation_row[GRN_inter_order,]
####oreder of whole GRS from best fit to worst fit
data_fit_order<-data_fit_order1[rep(1:9,159)+rep((order(apply(LL_list,1,min))-1)*9,each=9),]
annotation_row_order<-annotation_row_order1[rep(1:9,159)+rep((order(apply(LL_list,1,min))-1)*9,each=9),]
################################################################################
#############################Determine fitted genes#####################################
################################################################################
LL_logic_gates_threshold<-LL_logic_gates
LL_logic_gates_threshold[LL_logic_gates_threshold>fit_threshold]<-NA

LL_logic_gates_threshold[!is.na(LL_logic_gates_threshold)]<-1
LL_logic_gates_threshold[is.na(LL_logic_gates_threshold)]<-0
################################################################################
cols<- colorRampPalette(c("grey","gold2"))(2)

pdf("Exp_fitting_logics_threshold_order_plot.pdf",width=3,height=3)
pheatmap(LL_logic_gates_threshold[rev(do.call(order, as.data.frame(LL_logic_gates_threshold))),],
         cellheight=1, cellwidth = 6,gaps_col = c(3,9),
         legend=T,color = cols,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()
#############################count types of GRS###################################################
GRS_classification<-LL_logic_gates_threshold[rev(do.call(order, as.data.frame(LL_logic_gates_threshold))),]
GRS_Sum<-cbind(apply(GRS_classification[,1:3],1,sum),apply(GRS_classification[,4:9],1,sum),apply(GRS_classification[,10:17],1,sum))

single_TF_per<-length(which(GRS_Sum[,1]>0))/159
two_TF_per<-length(which(GRS_Sum[,2]>0 & GRS_Sum[,1]==0))/159
not_fit_per<-length(which(apply(GRS_Sum,1,sum)==0))/159
three_TF_per<-1-single_TF_per-two_TF_per-not_fit_per
#############################count types of GRS###################################################
Gate_count<-GRS_classification[which(apply(GRS_classification[,4:9],1,sum)>0),4:9]
Gate_count_total<-rbind(Gate_count[,c(1,2)],Gate_count[,c(3,4)],Gate_count[,c(5,6)])

AND_count<-length(intersect(which(Gate_count_total[,1]==0),which(Gate_count_total[,2]==1)))
OR_count<-length(intersect(which(Gate_count_total[,1]==1),which(Gate_count_total[,2]==0)))
ANDOR_count<-length(intersect(which(Gate_count_total[,1]==1),which(Gate_count_total[,2]==1)))

df<-data.frame(Gate_type=c("AND","OR","Both"),Count=c(AND_count,OR_count,ANDOR_count))
df$Gate_type <- factor(df$Gate_type,levels = c("AND","OR","Both"))

pdf("Gate_barplot.pdf",width=4,height=4)
ggplot(data=df, aes(x=Gate_type, y=Count)) +
  geom_bar(stat="identity",color="black")+
  geom_text(aes(label=Count), vjust=1.6, color="white", size=3.5)+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_minimal()
dev.off()

############################################################
############################################################
#########NFkB IRFs synergistically regulated genes##########
############################################################
############################################################
NFkBandIRF3_gene<-which(LL_logic_gates[,4]<fit_threshold)
tripleAND_gene<-which(LL_logic_gates[,10]<fit_threshold)
NFkBandIRF3orMAPK_gene<-which(LL_logic_gates[,13]<fit_threshold)
IRF3andNFkBorMAPK_gene<-which(LL_logic_gates[,15]<fit_threshold)
MAPKorIRF3andNFkB_gene<-which(LL_logic_gates[,16]<fit_threshold)

IRF3andNFkB_gene<-unique(c(names(NFkBandIRF3_gene),names(tripleAND_gene),names(NFkBandIRF3orMAPK_gene),
                           names(IRF3andNFkBorMAPK_gene),names(MAPKorIRF3andNFkB_gene)))

IRF3andNFkB_non_comp_gene<-unique(c(names(NFkBandIRF3_gene),names(tripleAND_gene)))

NFkBIRF3_reported<-c("Ccl5","Cxcl10","Gbp5","Irg1","Ifnb1")

load("IRF3KO_exp.RData")
WT_peak<-apply(Exp_data_mean[,1:4],1,max)
IRF3KO_peak<-apply(IRF3KO_exp,1,max)

IRF3_NFkB_nonComp_order<-pheatmap(log2(Exp_data_mean[IRF3andNFkB_non_comp_gene,c(1:8,13:16)]+0.01),cluster_cols = F,cluster_rows = T)
IRF3_NFkB_Comp_order<-pheatmap(log2(Exp_data_mean[setdiff(IRF3andNFkB_gene,IRF3andNFkB_non_comp_gene),c(1:8,13:16)]+0.01),cluster_cols = F,cluster_rows = T)

exp_genes<-c(which(names(WT_peak)==c("Ppp1r15a")),which(names(WT_peak)==c("Srgn")))
  
pdf("IRF3andNFkB_genes_all_display.pdf",width=4,height=4)
plot(log2(WT_peak+0.01),log2(IRF3KO_peak+0.01),pch=16,xlim=c(-7,9),ylim=c(-7,9),col="grey")
lines(log2(WT_peak[IRF3andNFkB_gene]+0.01),log2(IRF3KO_peak[IRF3andNFkB_gene]+0.01),pch=16,type="p",col="red")
lines(log2(WT_peak[NFkBIRF3_reported]+0.01),log2(IRF3KO_peak[NFkBIRF3_reported]+0.01),pch=1,type="p",col="black")
lines(log2(WT_peak[exp_genes]+0.01),log2(IRF3KO_peak[exp_genes]+0.01),pch=1,type="p",col="grey")
abline(coef = c(0,1), lty = 2,col="grey")
legend("topleft",legend=c("Identified","Previously reported"),
       col=c("red","black"),pch=c(16,1), cex=0.8)
graphics.off()
########################################################################################################################
########################################visualize for two individual genes##############################################
########################################################################################################################
RS_raw_plot<-RS_raw
RS_raw_plot[RS_raw_plot>=2]<-2
RS_raw_plot[RS_raw_plot<=(-1)]<-(-1)

annotation_row_fit_plot<-as.data.frame(cbind(logic_disp,RS_raw_plot)[,c(6,5,4,3,2,1)])
annotation_row_fit_plot<-cbind(target_fun_list_thre,annotation_row_fit_plot)
annotation_row_fit_plot[,2:4]<-RS_raw_plot[,c(3,2,1)]

annotation_row_plot<-as.data.frame(rbind(c(NA,NA,NA,NA,NA,NA,NA),c(NA,NA,NA,NA,NA,NA,NA))[rep(1,9*159),])
annotation_row_plot[c(1:dim(annotation_row_plot)[1])[-seq(1,dim(annotation_row_plot)[1],9)],]<-annotation_row_fit_plot

colnames(annotation_row_plot)<-c("LL","TF3","TF2","TF1","TF23","TF13","TF12")
rownames(annotation_row_plot)<-paste(rep(rownames(Exp_data_mean),each=9),rep(c("",c("  1a2a3","  1o2o3","  1o(2a3)","  1a(2o3)","  2o(1a3)","  2a(1o3)","  3o(1a2)","  3a(1o2)")),159),sep="")

anno_colors <- list("TF3"=c("firebrick","white"),
                    "TF2"=c("firebrick","white"),
                    "TF1"=c("firebrick","white"),
                    "TF23"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF13"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF12"=c("A"="#E69F00","O"="#56B4E9"),
                    "LL"=c("gold2", "grey"))

annotation_row_order1_plot<-annotation_row_plot[GRN_inter_order,]
annotation_row_order_plot<-annotation_row_order1_plot[rep(1:9,159)+rep((order(apply(LL_list,1,min))-1)*9,each=9),]
gene_ind<-c(which(rownames(data_fit_order)=="Srgn"),which(rownames(data_fit_order)=="Ppp1r15a"))

data_plot<-data_fit_order[rep(gene_ind,each=9)+seq(0,8,1),]
data_plot<-rbind(data_plot,"min"=data_plot[1,],"max"=data_plot[1,])

anno_plot<-annotation_row_order_plot[rep(gene_ind,each=9)+seq(0,8,1),]
anno_plot[which(anno_plot[,1]<fit_threshold),1]<-0
anno_plot[which(anno_plot[,1]>fit_threshold),1]<-1
range_keep<-as.data.frame(rbind("min"=c(annotation_row_order_plot[2,1],c(-1,-1,-1),c("A","A","A")),"max"=c(annotation_row_order_plot[2,1],c(2,2,2),c("A","A","A"))))
range_keep[,1:4]<-rbind(c(anno_plot[2,1],c(-1,-1,-1)),c(anno_plot[2,1],c(2,2,2)))
colnames(range_keep)<-colnames(anno_plot)
anno_plot<-rbind(anno_plot,range_keep)

pdf("Exp_fitting_two_genes_Srgn_Ppp1r15a.pdf",width=6,height=13)
pheatmap(log2(data_plot+0.01),cellheight=8, cellwidth = 3,gaps_row = c(9,18),gaps_col = seq(4,16,4),
         legend=T,
         show_rownames = T,show_colnames = F,cluster_cols = F,cluster_rows = F,font_size=0.1,
         annotation_row = anno_plot,annotation_legend = T, annotation_names_row = F,annotation_colors = anno_colors)
graphics.off()

data_plot_org<-data_plot
anno_plot_org<-anno_plot
anno_colors_org<-anno_colors
save(data_plot_org,anno_plot_org,anno_colors_org,file = "Two_genes_examples.RData")
