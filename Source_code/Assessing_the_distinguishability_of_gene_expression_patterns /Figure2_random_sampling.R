# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics")
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
#####################################################NFkB or(IRF3 and cJun)
#######################generate input
caRNA_time<-c(0,15,30,60)
Basal_level<-0
Medium_level<-0.5
input_raw<-as.data.frame(expand.grid(c(1,Medium_level,Basal_level),c(1,Medium_level,Basal_level),c(1,Medium_level,Basal_level)))[-27,]

SSE<-function(x,y,dynamics){
  return(sum((dynamics[x,]-dynamics[y,])^2))
}

sample_size<-1000
set.seed(18)
K1p<-runif(8*sample_size, min=-2, max=2)
K2p<-runif(8*sample_size, min=-2, max=2)
K3p<-runif(8*sample_size, min=-2, max=2)
K4p<-runif(8*sample_size, min=0, max=1)
K5p<-runif(8*sample_size, min=(-2), max=1)

K1p[sample(1:(8*sample_size),(8*sample_size)/10,replace=F)]<-6
K2p[sample(1:(8*sample_size),(8*sample_size)/10,replace=F)]<-6
K3p[sample(1:(8*sample_size),(8*sample_size)/10,replace=F)]<-6

RS_matrix<-rbind(c("None","None","None"),c("None","None","None"))[rep(1,8000),]
RS_matrix[which(10^(K1p)<=0.5),1]<-"S"
RS_matrix[intersect(which(10^(K1p)<=5),which(10^(K1p)>0.5)),1]<-"M"
RS_matrix[which(10^(K1p)>5),1]<-"W"
RS_matrix[which(K1p==6),1]<-"N"

RS_matrix[which(10^(K2p)<=0.5),2]<-"S"
RS_matrix[intersect(which(10^(K2p)<=5),which(10^(K2p)>0.5)),2]<-"M"
RS_matrix[which(10^(K2p)>5),2]<-"W"
RS_matrix[which(K2p==6),2]<-"N"

RS_matrix[which(10^(K3p)<=0.5),3]<-"S"
RS_matrix[intersect(which(10^(K3p)<=5),which(10^(K3p)>0.5)),3]<-"M"
RS_matrix[which(10^(K3p)>5),3]<-"W"
RS_matrix[which(K3p==6),3]<-"N"
##############################
############genereate data
##############################
dynamic_exp<-vector()
for(GRN_ind in 1:8){
  for(para_ind in 1:sample_size){
    Kd1<-10^(K1p[para_ind+sample_size*(GRN_ind-1)])
    Kd2<-10^(K2p[para_ind+sample_size*(GRN_ind-1)])
    Kd3<-10^(K3p[para_ind+sample_size*(GRN_ind-1)])
    ksyn<-1
    k0<-K4p[para_ind+sample_size*(GRN_ind-1)]
    kdeg<-10^(K5p[para_ind+sample_size*(GRN_ind-1)])
      
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
####################Normalize expression
dynamic_norm<-dynamic_exp*(100/apply(dynamic_exp[,1:4],1,max))
ksyn_list<-100/apply(dynamic_exp[,1:4],1,max)
induction<-apply(dynamic_norm[,1:4]/dynamic_norm[,1],1,function(x) any(x>=2))
dynamic_indcued<-dynamic_norm[induction,]
rownames(dynamic_indcued)<-c(1:dim(dynamic_indcued)[1])
########################################################################
########################################################################
########################################################################
logic_disp<-rbind(c("A","A","A"),c("O","O","O"),c("O","O","A"),c("A","A","O"),c("O","A","O"),c("A","O","A"),c("A","O","O"),c("O","A","A"))
annotation_row_raw<-cbind(logic_disp[rep(c(1:8),each=1000),],as.data.frame(RS_matrix,stringsAsFactors = T),K4p,log10(ksyn_list),K5p)[,c(9,8,7,6,5,4,3,2,1)]
annotation_row<-annotation_row_raw[induction,]

colfunc <- colorRampPalette(c("white", "firebrick"))
rownames(annotation_row) <- rownames(dynamic_indcued)
colnames(annotation_row)<-c("kdeg","ksyn","k0","TF3","TF2","TF1","TF23","TF13","TF12")
anno_colors <- list("kdeg"=c("forestgreen","white"),
                    "ksyn"= c("white", "orange"),
                    "k0"= c("white", "deepskyblue3"),
                    "TF3"=c("S"=colfunc(5)[5],"M"=colfunc(5)[3],"W"=colfunc(5)[2],"N"=colfunc(5)[1]),
                    "TF2"=c("S"=colfunc(5)[5],"M"=colfunc(5)[3],"W"=colfunc(5)[2],"N"=colfunc(5)[1]),
                    "TF1"=c("S"=colfunc(5)[5],"M"=colfunc(5)[3],"W"=colfunc(5)[2],"N"=colfunc(5)[1]),
                    "TF23"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF13"=c("A"="#E69F00","O"="#56B4E9"),
                    "TF12"=c("A"="#E69F00","O"="#56B4E9"))

con_reorder<-c(25,21,9,19,7,3,1)
con_ind_f<-unlist(lapply(con_reorder,function(x) seq((x*4-3),x*4,1)))

pdf("GRS_sampling.pdf",width=7,height=15)
pheatmap(dynamic_indcued[,con_ind_f],cluster_rows =T,cluster_cols = F,cellheight = 0.3,cellweight = 0.05,gaps_col =seq(4,24,4), 
         show_rownames = F,show_colnames = F,legend=T,clustering_method="single",
         annotation_row = annotation_row,annotation_legend = T, annotation_names_row = F,annotation_colors = anno_colors)
graphics.off()
########################################################################
SSE_7sets_sampling<-matrix(nrow=dim(dynamic_indcued)[1],ncol=dim(dynamic_indcued)[1])
dynamic_7sets<-dynamic_indcued[,con_ind_f]
for(i in 1:dim(dynamic_indcued)[1]){
  for(j in 1:dim(dynamic_indcued)[1]){
    SSE_7sets_sampling[i,j]<-SSE(i,j,dynamic_7sets)
  }
}
########################################################################
SSE_Amplitude_sampling<-matrix(nrow=dim(dynamic_indcued)[1],ncol=dim(dynamic_indcued)[1])
for(i in 1:dim(dynamic_indcued)[1]){
  for(j in 1:dim(dynamic_indcued)[1]){
    SSE_Amplitude_sampling[i,j]<-SSE(i,j,dynamic_indcued)
  }
}
################################################################################################
################################################################################################
########Validation for different dynamics###############################################################################################
################################################################################################
################################################################################################
source("Dyn_fit_eqa_k0.R")

SSE<-function(x,y,dynamics){
  return(sum((dynamics[x,]-dynamics[y,])^2))
}
############################################################################################
##############################################Delayed##############################################
############################################################################################
TF_time_ca<-c(0,10,11,30,31,40,41,60)

NFkB_H_org<-c(1,1,1,1,1,1,1,1)
NFkB_L_org<-c(Basal_level,Basal_level,Basal_level,Basal_level,1,1,1,1)
NFkB_N_org<-rep(Basal_level,8)

IRF3_H_org<-c(1,1,1,1,1,1,1,1)
IRF3_L_org<-c(Basal_level,Basal_level,Basal_level,Basal_level,1,1,1,1)
IRF3_N_org<-rep(Basal_level,8)

cJun_H_org<-c(1,1,1,1,1,1,1,1)
cJun_L_org<-c(Basal_level,Basal_level,Basal_level,Basal_level,1,1,1,1)
cJun_N_org<-rep(Basal_level,8)
######################################
int_length<-0.1
NFkB_H<-pchip(TF_time_ca, NFkB_H_org, seq(0,70,int_length))
IRF3_H<-pchip(TF_time_ca, IRF3_H_org, seq(0,70,int_length))
cJun_H<-pchip(TF_time_ca, cJun_H_org, seq(0,70,int_length))

NFkB_L<-pchip(TF_time_ca, NFkB_L_org, seq(0,70,int_length))
IRF3_L<-pchip(TF_time_ca, IRF3_L_org, seq(0,70,int_length))
cJun_L<-pchip(TF_time_ca, cJun_L_org, seq(0,70,int_length))

NFkB_N<-rep(Basal_level,(70/int_length+1))
IRF3_N<-rep(Basal_level,(70/int_length+1))
cJun_N<-rep(Basal_level,(70/int_length+1))
############################################################################
NFkB <- cbind(NFkB_H,NFkB_L,NFkB_N)[,rep(c(1,2,3),9)[-27]]
IRF3 <- cbind(IRF3_H,IRF3_L,IRF3_N)[,rep(rep(c(1,2,3),each=3),3)[-27]]
cJun <- cbind(cJun_H,cJun_L,cJun_N)[,rep(c(1,2,3),each=9)[-27]]
############################################################################
dynamic_exp<-vector()
for(GRN_ind in 1:8){
  for(para_ind in 1:sample_size){
    
    Kd1<-10^(K1p[para_ind+sample_size*(GRN_ind-1)])
    Kd2<-10^(K2p[para_ind+sample_size*(GRN_ind-1)])
    Kd3<-10^(K3p[para_ind+sample_size*(GRN_ind-1)])
    ksyn<-1
    k0<-K4p[para_ind+sample_size*(GRN_ind-1)]
    kdeg<-10^(K5p[para_ind+sample_size*(GRN_ind-1)])
    
    R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg
    
    Gene_exp<-matrix(ode(y =c(RNA1=R_basal,RNA2=R_basal,RNA3=R_basal,RNA4=R_basal,RNA5=R_basal,RNA6=R_basal,RNA7=R_basal,RNA8=R_basal,RNA9=R_basal,RNA10=R_basal,RNA11=R_basal,RNA12=R_basal,RNA13=R_basal,RNA14=R_basal,RNA15=R_basal,RNA16=R_basal,RNA17=R_basal,RNA18=R_basal,RNA19=R_basal,RNA20=R_basal,RNA21=R_basal,RNA22=R_basal,RNA23=R_basal,RNA24=R_basal,RNA25=R_basal,RNA26=R_basal), 
                         times = c(0,15,30,60),func =get(paste("P1_GRN",GRN_ind,"_fit",sep="")),parms = c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,ksyn=ksyn,kdeg=kdeg,k0=k0),hmax=1)[,-1],nrow=1)
    
    dynamic_exp<-rbind(dynamic_exp,Gene_exp)
  }
  print(GRN_ind)
}
dynamic_norm<-dynamic_exp*(100/apply(dynamic_exp[,1:4],1,max))
ksyn_list<-100/apply(dynamic_exp[,1:4],1,max)
dynamic_indcued<-dynamic_norm[which(induction),]
rownames(dynamic_indcued)<-c(1:dim(dynamic_indcued)[1])
#######################################################
SSE_all<-matrix(nrow=dim(dynamic_indcued)[1],ncol=dim(dynamic_indcued)[1])
for(i in 1:dim(dynamic_indcued)[1]){
  for(j in 1:dim(dynamic_indcued)[1]){
    SSE_all[i,j]<-SSE(i,j,dynamic_indcued)
  }
}

Delayed_dynamic<-dynamic_indcued
Delayed_SSE<-SSE_all
############################################################################################
##############################################Gradient##############################################
############################################################################################
int_length<-0.1
NFkB_L<-(1/60)*seq(0,70,int_length)
IRF3_L<-(1/60)*seq(0,70,int_length)
cJun_L<-(1/60)*seq(0,70,int_length)
############################################################################
NFkB <- cbind(NFkB_H,NFkB_L,NFkB_N)[,rep(c(1,2,3),9)[-27]]
IRF3 <- cbind(IRF3_H,IRF3_L,IRF3_N)[,rep(rep(c(1,2,3),each=3),3)[-27]]
cJun <- cbind(cJun_H,cJun_L,cJun_N)[,rep(c(1,2,3),each=9)[-27]]

dynamic_exp<-vector()
for(GRN_ind in 1:8){
  for(para_ind in 1:sample_size){
    
    Kd1<-10^(K1p[para_ind+sample_size*(GRN_ind-1)])
    Kd2<-10^(K2p[para_ind+sample_size*(GRN_ind-1)])
    Kd3<-10^(K3p[para_ind+sample_size*(GRN_ind-1)])
    ksyn<-1
    k0<-K4p[para_ind+sample_size*(GRN_ind-1)]
    kdeg<-10^(K5p[para_ind+sample_size*(GRN_ind-1)])
    
    R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg
    
    Gene_exp<-matrix(ode(y =c(RNA1=R_basal,RNA2=R_basal,RNA3=R_basal,RNA4=R_basal,RNA5=R_basal,RNA6=R_basal,RNA7=R_basal,RNA8=R_basal,RNA9=R_basal,RNA10=R_basal,RNA11=R_basal,RNA12=R_basal,RNA13=R_basal,RNA14=R_basal,RNA15=R_basal,RNA16=R_basal,RNA17=R_basal,RNA18=R_basal,RNA19=R_basal,RNA20=R_basal,RNA21=R_basal,RNA22=R_basal,RNA23=R_basal,RNA24=R_basal,RNA25=R_basal,RNA26=R_basal), 
                         times = c(0,15,30,60),func =get(paste("P1_GRN",GRN_ind,"_fit",sep="")),parms = c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,ksyn=ksyn,kdeg=kdeg,k0=k0),hmax=1)[,-1],nrow=1)
    
    dynamic_exp<-rbind(dynamic_exp,Gene_exp)
  }
  print(GRN_ind)
}
dynamic_norm<-dynamic_exp*(100/apply(dynamic_exp[,1:4],1,max))
ksyn_list<-100/apply(dynamic_exp[,1:4],1,max)
dynamic_indcued<-dynamic_norm[which(induction),]
rownames(dynamic_indcued)<-c(1:dim(dynamic_indcued)[1])
#######################################################
SSE_all<-matrix(nrow=dim(dynamic_indcued)[1],ncol=dim(dynamic_indcued)[1])
for(i in 1:dim(dynamic_indcued)[1]){
  for(j in 1:dim(dynamic_indcued)[1]){
    SSE_all[i,j]<-SSE(i,j,dynamic_indcued)
  }
}
Gradient_dynamic<-dynamic_indcued
Gradient_SSE<-SSE_all
############################################################################################
##############################################Transient##############################################
############################################################################################
int_length<-0.1
NFkB_L_org<-c(1,1,1,1,Basal_level,Basal_level,Basal_level,Basal_level)
IRF3_L_org<-c(1,1,1,1,Basal_level,Basal_level,Basal_level,Basal_level)
cJun_L_org<-c(1,1,1,1,Basal_level,Basal_level,Basal_level,Basal_level)
####################################################################################
NFkB_L<-pchip(TF_time_ca, NFkB_L_org, seq(0,70,int_length))
IRF3_L<-pchip(TF_time_ca, IRF3_L_org, seq(0,70,int_length))
cJun_L<-pchip(TF_time_ca, cJun_L_org, seq(0,70,int_length))
############################################################################
NFkB <- cbind(NFkB_H,NFkB_L,NFkB_N)[,rep(c(1,2,3),9)[-27]]
IRF3 <- cbind(IRF3_H,IRF3_L,IRF3_N)[,rep(rep(c(1,2,3),each=3),3)[-27]]
cJun <- cbind(cJun_H,cJun_L,cJun_N)[,rep(c(1,2,3),each=9)[-27]]

NFkB_list<-cbind("H","L","N")[,rep(c(1,2,3),9)[-27]]
IRF3_list<-cbind("H","L","N")[,rep(rep(c(1,2,3),each=3),3)[-27]]
cJun_list <- cbind("H","L","N")[,rep(c(1,2,3),each=9)[-27]]
############################################################################
dynamic_exp<-vector()
for(GRN_ind in 1:8){
  for(para_ind in 1:sample_size){
    
    Kd1<-10^(K1p[para_ind+sample_size*(GRN_ind-1)])
    Kd2<-10^(K2p[para_ind+sample_size*(GRN_ind-1)])
    Kd3<-10^(K3p[para_ind+sample_size*(GRN_ind-1)])
    ksyn<-1
    k0<-K4p[para_ind+sample_size*(GRN_ind-1)]
    kdeg<-10^(K5p[para_ind+sample_size*(GRN_ind-1)])
    
    R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg
    
    Gene_exp<-matrix(ode(y =c(RNA1=R_basal,RNA2=R_basal,RNA3=R_basal,RNA4=R_basal,RNA5=R_basal,RNA6=R_basal,RNA7=R_basal,RNA8=R_basal,RNA9=R_basal,RNA10=R_basal,RNA11=R_basal,RNA12=R_basal,RNA13=R_basal,RNA14=R_basal,RNA15=R_basal,RNA16=R_basal,RNA17=R_basal,RNA18=R_basal,RNA19=R_basal,RNA20=R_basal,RNA21=R_basal,RNA22=R_basal,RNA23=R_basal,RNA24=R_basal,RNA25=R_basal,RNA26=R_basal), 
                         times = c(0,15,30,60),func =get(paste("P1_GRN",GRN_ind,"_fit",sep="")),parms = c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,ksyn=ksyn,kdeg=kdeg,k0=k0),hmax=1)[,-1],nrow=1)
    
    dynamic_exp<-rbind(dynamic_exp,Gene_exp)
  }
  print(GRN_ind)
}
dynamic_norm<-dynamic_exp*(100/apply(dynamic_exp[,1:4],1,max))
ksyn_list<-100/apply(dynamic_exp[,1:4],1,max)
dynamic_indcued<-dynamic_norm[which(induction),]
rownames(dynamic_indcued)<-c(1:dim(dynamic_indcued)[1])
#######################################################
SSE_all<-matrix(nrow=dim(dynamic_indcued)[1],ncol=dim(dynamic_indcued)[1])
for(i in 1:dim(dynamic_indcued)[1]){
  for(j in 1:dim(dynamic_indcued)[1]){
    SSE_all[i,j]<-SSE(i,j,dynamic_indcued)
  }
}

Transient_dynamic<-dynamic_indcued
Transient_SSE<-SSE_all
############################################################################
SSE_7sets_sampling[SSE_7sets_sampling<0.1]<-0
hclust7_single <-hclust(as.dist(SSE_7sets_sampling), method = "single", members = NULL)
dend7 <- hclust7_single %>% as.dendrogram
##########################################################################
##########################################################################
##########################################################################
SSE_Amplitude_sampling[SSE_Amplitude_sampling<0.1]<-0

hclustall_single<-hclust(as.dist(SSE_Amplitude_sampling), method = "single", members = NULL)
dendall <-hclustall_single%>%as.dendrogram
hclust_27HML<-hclustall_single
##########################################################################
Delayed_SSE[Delayed_SSE<0.1]<-0

hclustdelayed_single<-hclust(as.dist(Delayed_SSE), method = "single", members = NULL)
dend_delayed<-hclustdelayed_single%>%as.dendrogram
##########################################################################
Gradient_SSE[Gradient_SSE<0.1]<-0

hclustgradient_single<-hclust(as.dist(Gradient_SSE), method = "single", members = NULL)
dend_gradient<-hclustgradient_single%>%as.dendrogram
##########################################################################
Transient_SSE[Transient_SSE<0.1]<-0

hclusttransient_single<-hclust(as.dist(Transient_SSE), method = "single", members = NULL)
dend_transient<-hclusttransient_single%>%as.dendrogram
##########################################################################
##########################################################################threshold-cluster number plot
##########################################################################
duplicated_id_7<-which(duplicated(names(heights_per_k.dendrogram(dend7))))
name_list_7<-heights_per_k.dendrogram(dend7)
name_list_nonR_7<-name_list_7
if(length(duplicated_id_7)!=0){
  name_list_nonR_7<-name_list_7[-duplicated_id_7]
}

if(any(name_list_nonR_7<0.1)){
  name_list_nonR_7<-name_list_nonR_7[-which(name_list_nonR_7<0.1)]
}

den_plot_7<-name_list_nonR_7
den_plot_select_log10_7<-round(log10(den_plot_7),digits = 2)[-1]

plot_thre_7<-cbind(seq(-2,5,0.01),rep(0,length(seq(-2,5,0.01))))

for(i in 1:length(den_plot_select_log10_7)){
  ind_s<-intersect(which(plot_thre_7[,1]<=den_plot_select_log10_7[i]),which(plot_thre_7[,1]>den_plot_select_log10_7[i+1]))
  plot_thre_7[ind_s,2]<-as.numeric(names(den_plot_select_log10_7)[i])
}
plot_thre_7[which(plot_thre_7[,1]>=den_plot_select_log10_7[1]),2]<-as.numeric(names(den_plot_select_log10_7)[1])
plot_thre_7[which(plot_thre_7[,1]<=tail(den_plot_select_log10_7,n=1)),2]<-as.numeric(tail(names(den_plot_select_log10_7),n=1))
##########################################################################
##########################################################################
##########################################################################All 26
##########################################################################
duplicated_id_all<-which(duplicated(names(heights_per_k.dendrogram(dendall))))
name_list_all<-heights_per_k.dendrogram(dendall)
name_list_nonR_all<-name_list_all
if(length(duplicated_id_all)!=0){
  name_list_nonR_all<-name_list_all[-duplicated_id_all]
}
if(any(name_list_nonR_all<0.1)){
  name_list_nonR_all<-name_list_nonR_all[-which(name_list_nonR_all<0.1)]
}

den_plot_all<-name_list_nonR_all
den_plot_select_log10_all<-round(log10(den_plot_all),digits = 2)[-1]

plot_thre_all<-cbind(seq(-2,5,0.01),rep(0,length(seq(-2,5,0.01))))

for(i in 1:length(den_plot_select_log10_all)){
  ind_s<-intersect(which(plot_thre_all[,1]<=den_plot_select_log10_all[i]),which(plot_thre_all[,1]>den_plot_select_log10_all[i+1]))
  plot_thre_all[ind_s,2]<-as.numeric(names(den_plot_select_log10_all)[i])
}
plot_thre_all[which(plot_thre_all[,1]>=den_plot_select_log10_all[1]),2]<-as.numeric(names(den_plot_select_log10_all)[1])
plot_thre_all[which(plot_thre_all[,1]<=tail(den_plot_select_log10_all,n=1)),2]<-as.numeric(tail(names(den_plot_select_log10_all),n=1))
##########################################################################
##########################################################################All delayed
duplicated_id_delayed<-which(duplicated(names(heights_per_k.dendrogram(dend_delayed))))
name_list_delayed<-heights_per_k.dendrogram(dend_delayed)
name_list_nonR_delayed<-name_list_delayed
if(length(duplicated_id_delayed)!=0){
  name_list_nonR_delayed<-name_list_delayed[-duplicated_id_delayed]
}
if(any(name_list_nonR_delayed<0.1)){
  name_list_nonR_delayed<-name_list_nonR_delayed[-which(name_list_nonR_delayed<0.1)]
}

den_plot_delayed<-name_list_nonR_delayed
den_plot_select_log10_delayed<-round(log10(den_plot_delayed),digits = 2)[-1]

plot_thre_delayed<-cbind(seq(-2,5,0.01),rep(0,length(seq(-2,5,0.01))))

for(i in 1:length(den_plot_select_log10_delayed)){
  ind_s<-intersect(which(plot_thre_delayed[,1]<=den_plot_select_log10_delayed[i]),which(plot_thre_delayed[,1]>den_plot_select_log10_delayed[i+1]))
  plot_thre_delayed[ind_s,2]<-as.numeric(names(den_plot_select_log10_delayed)[i])
}
plot_thre_delayed[which(plot_thre_delayed[,1]>=den_plot_select_log10_delayed[1]),2]<-as.numeric(names(den_plot_select_log10_delayed)[1])
plot_thre_delayed[which(plot_thre_delayed[,1]<=tail(den_plot_select_log10_delayed,n=1)),2]<-as.numeric(tail(names(den_plot_select_log10_delayed),n=1))
##########################################################################
##########################################################################All gradient
duplicated_id_gradient<-which(duplicated(names(heights_per_k.dendrogram(dend_gradient))))
name_list_gradient<-heights_per_k.dendrogram(dend_gradient)
name_list_nonR_gradient<-name_list_gradient
if(length(duplicated_id_gradient)!=0){
  name_list_nonR_gradient<-name_list_gradient[-duplicated_id_gradient]
}
if(any(name_list_nonR_gradient<0.1)){
  name_list_nonR_gradient<-name_list_nonR_gradient[-which(name_list_nonR_gradient<0.1)]
}

den_plot_gradient<-name_list_nonR_gradient
den_plot_select_log10_gradient<-round(log10(den_plot_gradient),digits = 2)[-1]

plot_thre_gradient<-cbind(seq(-2,5,0.01),rep(0,length(seq(-2,5,0.01))))

for(i in 1:length(den_plot_select_log10_gradient)){
  ind_s<-intersect(which(plot_thre_gradient[,1]<=den_plot_select_log10_gradient[i]),which(plot_thre_gradient[,1]>den_plot_select_log10_gradient[i+1]))
  plot_thre_gradient[ind_s,2]<-as.numeric(names(den_plot_select_log10_gradient)[i])
}
plot_thre_gradient[which(plot_thre_gradient[,1]>=den_plot_select_log10_gradient[1]),2]<-as.numeric(names(den_plot_select_log10_gradient)[1])
plot_thre_gradient[which(plot_thre_gradient[,1]<=tail(den_plot_select_log10_gradient,n=1)),2]<-as.numeric(tail(names(den_plot_select_log10_gradient),n=1))
############################################################
####################All transient###########################
duplicated_id_transient<-which(duplicated(names(heights_per_k.dendrogram(dend_transient))))
name_list_transient<-heights_per_k.dendrogram(dend_transient)
name_list_nonR_transient<-name_list_transient
if(length(duplicated_id_transient)!=0){
  name_list_nonR_transient<-name_list_transient[-duplicated_id_transient]
}
if(any(name_list_nonR_transient<0.1)){
  name_list_nonR_transient<-name_list_nonR_transient[-which(name_list_nonR_transient<0.1)]
}

den_plot_transient<-name_list_nonR_transient
den_plot_select_log10_transient<-round(log10(den_plot_transient),digits = 2)[-1]

plot_thre_transient<-cbind(seq(-2,5,0.01),rep(0,length(seq(-2,5,0.01))))

for(i in 1:length(den_plot_select_log10_transient)){
  ind_s<-intersect(which(plot_thre_transient[,1]<=den_plot_select_log10_transient[i]),which(plot_thre_transient[,1]>den_plot_select_log10_transient[i+1]))
  plot_thre_transient[ind_s,2]<-as.numeric(names(den_plot_select_log10_transient)[i])
}
plot_thre_transient[which(plot_thre_transient[,1]>=den_plot_select_log10_transient[1]),2]<-as.numeric(names(den_plot_select_log10_transient)[1])
plot_thre_transient[which(plot_thre_transient[,1]<=tail(den_plot_select_log10_transient,n=1)),2]<-as.numeric(tail(names(den_plot_select_log10_transient),n=1))
##########################

pdf("Cluster_treeheight_single_MDTG_sampling.pdf",width=6,height=6)
plot(plot_thre_7[101:701,1],plot_thre_7[101:701,2],col="black",xlab="",ylab="",type="o",cex=1.5,cex.axis=1.5,pch=16,lwd=3)
lines(plot_thre_all[101:701,1],plot_thre_all[101:701,2],col=brewer.pal(4, "Set1")[1],type="o",pch=16,cex=1.5,lwd=3)
lines(plot_thre_delayed[101:701,1],plot_thre_delayed[101:701,2],col=brewer.pal(4, "Set1")[2],type="o",pch=16,cex=1.5,lwd=3)
lines(plot_thre_gradient[101:701,1],plot_thre_gradient[101:701,2],col=brewer.pal(4, "Set1")[3],type="o",pch=16,cex=1.5,lwd=3)
lines(plot_thre_transient[101:701,1],plot_thre_transient[101:701,2],col=brewer.pal(4, "Set1")[4],type="o",pch=16,cex=1.5,lwd=3)
# legend("topright",legend=c("7 Datasets","Amplitute","Delayed","Gradient","Transient"),
#        col=c("black",brewer.pal(4, "Set1")),pch=c(16,16,16,16,16), cex=1)
graphics.off()
