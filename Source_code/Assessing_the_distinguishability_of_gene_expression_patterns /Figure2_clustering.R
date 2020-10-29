# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics/Submit_code/Source_code/Assessing_GRS_identifiability")
rm(list=ls())
library("MASS")
library(deSolve)
library(signal)
library(fields)
library(pheatmap)
library(RColorBrewer)
library(cluster)
library(reshape)
library(ggplot2)
require(cowplot)
source("Nondyn_opt_newbasalnew.R")
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################NFkB or(IRF3 and ISGF3)
#######################generate input
caRNA_time<-c(0,15,30,60)
Basal_level<-0
k0<-0
SSE<-function(x,y,dynamics){
  return(sum((dynamics[x,]-dynamics[y,])^2))
}
Medium_level<-0.5
input_raw<-as.data.frame(expand.grid(c(1,Medium_level,Basal_level),c(1,Medium_level,Basal_level),c(1,Medium_level,Basal_level)))[-27,]
para_list<-as.data.frame(expand.grid(c(0.1,1,10),c(0.1,1,10),c(0.1,1,10))[-c(27),])
##############################
############genereate data
##############################
dynamic_exp<-vector()
for(GRN_ind in 1:8){
  for(para_ind in 1:dim(para_list)[1]){
    Kd1<-para_list[para_ind,1]
    Kd2<-para_list[para_ind,2]
    Kd3<-para_list[para_ind,3]
    ksyn<-1
    kdeg<-0.1
    
    parS<-c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,k0=k0,ksyn=ksyn,kdeg=kdeg)
    R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg
    Gene_exp<-c(get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[1,1],h2=input_raw[1,2],h3=input_raw[1,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[2,1],h2=input_raw[2,2],h3=input_raw[2,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[3,1],h2=input_raw[3,2],h3=input_raw[3,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[4,1],h2=input_raw[4,2],h3=input_raw[4,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[5,1],h2=input_raw[5,2],h3=input_raw[5,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[6,1],h2=input_raw[6,2],h3=input_raw[6,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[7,1],h2=input_raw[7,2],h3=input_raw[7,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[8,1],h2=input_raw[8,2],h3=input_raw[8,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[9,1],h2=input_raw[9,2],h3=input_raw[9,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[10,1],h2=input_raw[10,2],h3=input_raw[10,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[11,1],h2=input_raw[11,2],h3=input_raw[11,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[12,1],h2=input_raw[12,2],h3=input_raw[12,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[13,1],h2=input_raw[13,2],h3=input_raw[13,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[14,1],h2=input_raw[14,2],h3=input_raw[14,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[15,1],h2=input_raw[15,2],h3=input_raw[15,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[16,1],h2=input_raw[16,2],h3=input_raw[16,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[17,1],h2=input_raw[17,2],h3=input_raw[17,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[18,1],h2=input_raw[18,2],h3=input_raw[18,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[19,1],h2=input_raw[19,2],h3=input_raw[19,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[20,1],h2=input_raw[20,2],h3=input_raw[20,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[21,1],h2=input_raw[21,2],h3=input_raw[21,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[22,1],h2=input_raw[22,2],h3=input_raw[22,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[23,1],h2=input_raw[23,2],h3=input_raw[23,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[24,1],h2=input_raw[24,2],h3=input_raw[24,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[25,1],h2=input_raw[25,2],h3=input_raw[25,3])),
                get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[26,1],h2=input_raw[26,2],h3=input_raw[26,3])))
    dynamic_exp<-rbind(dynamic_exp,Gene_exp)
  }
}
logic_name<-c("1a2a3","1o2o3","1o(2a3)","1a(2o3)","2o(1a3)","2a(1o3)","3o(1a2)","3a(1o2)")
para_name<-paste(expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,1],
                 expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,2],
                 expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,3],sep="")[-27]
GRS_name<-paste(rep(logic_name,each=26),rep(para_name,8))
rownames(dynamic_exp)<-GRS_name
########################
########################gate select for weak GRS
G_list<-vector()
for(GRN_ind in 1:8){
  for(para_ind in 1:26){
    Kd1<-para_list[para_ind,1]
    Kd2<-para_list[para_ind,2]
    Kd3<-para_list[para_ind,3]
    
    parS<-c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,k0=k0)
    
    G_exp<-c(get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[1,1],h2=input_raw[1,2],h3=input_raw[1,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[2,1],h2=input_raw[2,2],h3=input_raw[2,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[3,1],h2=input_raw[3,2],h3=input_raw[3,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[4,1],h2=input_raw[4,2],h3=input_raw[4,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[5,1],h2=input_raw[5,2],h3=input_raw[5,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[6,1],h2=input_raw[6,2],h3=input_raw[6,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[7,1],h2=input_raw[7,2],h3=input_raw[7,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[8,1],h2=input_raw[8,2],h3=input_raw[8,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[9,1],h2=input_raw[9,2],h3=input_raw[9,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[10,1],h2=input_raw[10,2],h3=input_raw[10,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[11,1],h2=input_raw[11,2],h3=input_raw[11,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[12,1],h2=input_raw[12,2],h3=input_raw[12,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[13,1],h2=input_raw[13,2],h3=input_raw[13,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[14,1],h2=input_raw[14,2],h3=input_raw[14,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[15,1],h2=input_raw[15,2],h3=input_raw[15,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[16,1],h2=input_raw[16,2],h3=input_raw[16,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[17,1],h2=input_raw[17,2],h3=input_raw[17,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[18,1],h2=input_raw[18,2],h3=input_raw[18,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[19,1],h2=input_raw[19,2],h3=input_raw[19,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[20,1],h2=input_raw[20,2],h3=input_raw[20,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[21,1],h2=input_raw[21,2],h3=input_raw[21,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[22,1],h2=input_raw[22,2],h3=input_raw[22,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[23,1],h2=input_raw[23,2],h3=input_raw[23,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[24,1],h2=input_raw[24,2],h3=input_raw[24,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[25,1],h2=input_raw[25,2],h3=input_raw[25,3])),
             get(paste("G_GRN",GRN_ind,sep=""))(c(parS,h1=input_raw[26,1],h2=input_raw[26,2],h3=input_raw[26,3])))
    G_list<-rbind(G_list,G_exp)
  }
}
##remove weak GRS, 1:61
max_G<-apply(G_list,1,max)
weak_list<-order(max_G)[1:61]
dynamic_active<-dynamic_exp[-weak_list,]
# save(weak_list,file="weak_list.RData")
################################################
################################################
nonR_1o<-c(40,41,48:55)
nonR_1a<-c(60:61,66:71)

nonR_2o<-c(74,77,82,85,88:93)
nonR_2a<-c(96,99,102,105:109)

nonR_3o<-c(112,115:118,121,124:127)
nonR_3a<-c(134,137:139,142,145:147)

dynamic_nonR<-dynamic_active[-c(nonR_1o,nonR_1a,nonR_2o,nonR_2a,nonR_3o,nonR_3a),]
####################Normalize expression
dynamic_norm<-dynamic_nonR*(100/apply(dynamic_nonR,1,max))
dynamic_norm_basalN<-dynamic_norm
################################################################################################
################################################################################################
#################################Calculate Identifiability for 7 sets###################
################################################################################################
################################################################################################
con_reorder<-c(25,21,9,19,7,3,1)
con_ind_f<-unlist(lapply(con_reorder,function(x) seq((x*4-3),x*4,1)))
gene_exp<-dynamic_norm[,con_ind_f]
con_ind_f<-unlist(lapply(con_reorder,function(x) seq((x*4-3),x*4,1)))
SSE_7sets<-matrix(nrow=dim(dynamic_norm)[1],ncol=dim(dynamic_norm)[1])
for(i in 1:dim(dynamic_norm)[1]){
  for(j in 1:dim(dynamic_norm)[1]){
    SSE_7sets[i,j]<-SSE(i,j,dynamic_norm[,con_ind_f])
  }
}
Gene_exp_Medium<-dynamic_norm_basalN
Identifiability_7sets_SSE_basalN<-apply(SSE_7sets,1,function(x) sort(x)[2])
# save(Gene_exp_Medium,SSE_7sets,Identifiability_7sets_SSE_basalN,dynamic_norm,file="Identifiability_7sets_SSE_basalN.RData")
# load("Identifiability_7sets_SSE_basalN.RData")
################################################################################################
#######################################Gene expression cluster##################################
################################################################################################
colfunc <- colorRampPalette(c("white", "firebrick"))
logic_disp<-rbind(c("A","A","A"),c("O","O","O"),c("O","O","A"),c("A","A","O"),c("O","A","O"),c("A","O","A"),c("A","O","O"),c("O","A","A"))

annotation_row<-cbind(logic_disp[rep(c(1:8),each=26),],expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[-27,][rep(c(1:26),8),])[rownames(dynamic_exp)%in%rownames(dynamic_norm),c(6,5,4,3,2,1)]

rownames(annotation_row) <- rownames(gene_exp)
colnames(annotation_row)<-c("TF3","TF2","TF1","TF23","TF13","TF12")
anno_colors <- list("TF3"=c("S"=colfunc(5)[5],"M"=colfunc(5)[3],"W"=colfunc(5)[2]),
                    "TF2"=c("S"=colfunc(5)[5],"M"=colfunc(5)[3],"W"=colfunc(5)[2]),
                    "TF1"=c("S"=colfunc(5)[5],"M"=colfunc(5)[3],"W"=colfunc(5)[2]),
                    "TF23"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF13"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF12"=c("A"="#E69F00","O"="#56B4E9"))
##############################################################
pdf("Nondyn_geneexp_clustering_basalN_single_plot.pdf",width=10,height=8)
pheatmap(gene_exp,cellheight=5, cellwidth = 10,legend=F,gaps_col = c(4,8,12,16,20,24),clustering_distance_rows =as.dist(SSE_7sets),
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = T,clustering_method="single",
         annotation_row = annotation_row,annotation_legend = T, annotation_names_row = F,annotation_colors = anno_colors)
graphics.off()

hclust_7sets_single<-pheatmap(gene_exp,cellheight=10, cellwidth = 10,legend=F,
                              show_rownames = T,show_colnames = F,cluster_cols = F,cluster_rows = T,clustering_distance_rows =as.dist(SSE_7sets),clustering_method="single",
                              annotation_row = annotation_row,annotation_legend = F, annotation_names_row = F,annotation_colors = anno_colors)
save(hclust_7sets_single,file="hclust_7sets_single.RData")
##############################################################
pdf("Nondyn_geneexp_clustering_basalN_complete.pdf",width=10,height=15)
pheatmap(gene_exp,cellheight=5, cellwidth = 10,legend=F,gaps_col = c(4,8,12,16,20,24),clustering_distance_rows =as.dist(SSE_7sets),         
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = T,clustering_method="complete",
         annotation_row = annotation_row,annotation_legend = F, annotation_names_row = F,annotation_colors = anno_colors)
graphics.off()
##############################################################
pdf("Nondyn_geneexp_clustering_basalN_average.pdf",width=10,height=15)
pheatmap(gene_exp,cellheight=5, cellwidth = 10,legend=F,gaps_col = c(4,8,12,16,20,24),clustering_distance_rows =as.dist(SSE_7sets),        
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = T,clustering_method="average",
         annotation_row = annotation_row,annotation_legend = F, annotation_names_row = F,annotation_colors = anno_colors)
graphics.off()
