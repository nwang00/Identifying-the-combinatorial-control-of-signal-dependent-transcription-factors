# setwd("/Users/ning/GoogleDrive/work/gene_regulation_network/code_dynamics")
rm(list=ls())
library(foreach)
library(iterators)
library(parallel)
library(plyr)
library(doMC)
library("MASS")
library(deSolve)
library(signal)
library(pheatmap)
library(RColorBrewer)
source("Nondyn_opt_newbasalnew.R")
doMC::registerDoMC(cores=20)
deltaN<-1.05
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################NFkB or(IRF3 and ISGF3)
#######################generate input
# caRNA_time<-c(0,0,4,13,15,17,28,30,32,56,60,60)
# caRNA_time<-c(0,15,30,60)
int_length<-0.1
caRNA_time<-seq(0,80,int_length)

Basal_level<-0
Medium_strength<-0.5

input_raw<-as.data.frame(expand.grid(c(1,Medium_strength,Basal_level),c(1,Medium_strength,Basal_level),c(1,Medium_strength,Basal_level)))[-27,]
para_list<-as.data.frame(expand.grid(c(0.1,1,10),c(0.1,1,10),c(0.1,1,10)))[-27,]
rep_size<-1000
######remove all weak
dynamic_exp_true<-vector()
dynamic_exp_rep<-vector()
for(GRN_ind in 1:8){
  for(para_ind in 1:dim(para_list)[1]){
    Kd1<-para_list[para_ind,1]
    Kd2<-para_list[para_ind,2]
    Kd3<-para_list[para_ind,3]
    ksyn<-1
    kdeg<-0.1
    k0<-0
    
    parS<-c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,k0=k0,ksyn=ksyn,kdeg=kdeg)
    R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg
    Gene_exp_true<-c(get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_raw[1,1],h2=input_raw[1,2],h3=input_raw[1,3])),
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
    dynamic_exp_true<-rbind(dynamic_exp_true,Gene_exp_true)
  }
}

set.seed(1)
dynamic_exp_rep<-list()
para_noise<-rnorm(5*rep_size*26*8*26, mean = 0, sd = deltaN*0.02)
k0_noise<-runif(rep_size*26*8*26, min = 0, max = deltaN*0.02)
input_error<-matrix(rnorm(26*3*rep_size*26*8,mean=0,sd=deltaN*0.02),26,3*rep_size*26*8) 
input_error[input_error<(-1)]<-0

dynamic_exp_rep<-foreach(rep_ind=1:rep_size) %dopar% {
  GRS_all<-vector()
  noise_ind<-(rep_ind-1)*208+1
  for(GRN_ind in 1:8){
    for(para_ind in 1:dim(para_list)[1]){
      #########adding noise
      Kd1<-para_list[para_ind,1]*exp(para_noise[(26*(5*noise_ind-4)-25):(26*(5*noise_ind-4))])
      Kd2<-para_list[para_ind,2]*exp(para_noise[(26*(5*noise_ind-3)-25):(26*(5*noise_ind-3))])
      Kd3<-para_list[para_ind,3]*exp(para_noise[(26*(5*noise_ind-2)-25):(26*(5*noise_ind-2))])
      
      ksyn<-1*exp(para_noise[(26*(5*noise_ind-1)-25):(26*(5*noise_ind-1))])
      kdeg<-0.1*exp(para_noise[(130*noise_ind-25):(130*noise_ind)])
      k0<-0+k0_noise[(26*noise_ind-25):(26*noise_ind)]
      
      input_noisy<-input_raw*(1+input_error[,(noise_ind*3-2):(noise_ind*3)])
      #####add 1
      noise_ind<-noise_ind+1
      
      R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg
      
      Gene_exp_rep<-c(get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[1],c(Kd1=Kd1[1],Kd2=Kd2[1],Kd3=Kd3[1],k0=k0[1],ksyn=ksyn[1],kdeg=kdeg[1],h1=input_noisy[1,1],h2=input_noisy[1,2],h3=input_noisy[1,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[2],c(Kd1=Kd1[2],Kd2=Kd2[2],Kd3=Kd3[2],k0=k0[2],ksyn=ksyn[2],kdeg=kdeg[2],h1=input_noisy[2,1],h2=input_noisy[2,2],h3=input_noisy[2,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[3],c(Kd1=Kd1[3],Kd2=Kd2[3],Kd3=Kd3[3],k0=k0[3],ksyn=ksyn[3],kdeg=kdeg[3],h1=input_noisy[3,1],h2=input_noisy[3,2],h3=input_noisy[3,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[4],c(Kd1=Kd1[4],Kd2=Kd2[4],Kd3=Kd3[4],k0=k0[4],ksyn=ksyn[4],kdeg=kdeg[4],h1=input_noisy[4,1],h2=input_noisy[4,2],h3=input_noisy[4,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[5],c(Kd1=Kd1[5],Kd2=Kd2[5],Kd3=Kd3[5],k0=k0[5],ksyn=ksyn[5],kdeg=kdeg[5],h1=input_noisy[5,1],h2=input_noisy[5,2],h3=input_noisy[5,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[6],c(Kd1=Kd1[6],Kd2=Kd2[6],Kd3=Kd3[6],k0=k0[6],ksyn=ksyn[6],kdeg=kdeg[6],h1=input_noisy[6,1],h2=input_noisy[6,2],h3=input_noisy[6,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[7],c(Kd1=Kd1[7],Kd2=Kd2[7],Kd3=Kd3[7],k0=k0[7],ksyn=ksyn[7],kdeg=kdeg[7],h1=input_noisy[7,1],h2=input_noisy[7,2],h3=input_noisy[7,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[8],c(Kd1=Kd1[8],Kd2=Kd2[8],Kd3=Kd3[8],k0=k0[8],ksyn=ksyn[8],kdeg=kdeg[8],h1=input_noisy[8,1],h2=input_noisy[8,2],h3=input_noisy[8,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[9],c(Kd1=Kd1[9],Kd2=Kd2[9],Kd3=Kd3[9],k0=k0[9],ksyn=ksyn[9],kdeg=kdeg[9],h1=input_noisy[9,1],h2=input_noisy[9,2],h3=input_noisy[9,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[10],c(Kd1=Kd1[10],Kd2=Kd2[10],Kd3=Kd3[10],k0=k0[10],ksyn=ksyn[10],kdeg=kdeg[10],h1=input_noisy[10,1],h2=input_noisy[10,2],h3=input_noisy[10,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[11],c(Kd1=Kd1[11],Kd2=Kd2[11],Kd3=Kd3[11],k0=k0[11],ksyn=ksyn[11],kdeg=kdeg[11],h1=input_noisy[11,1],h2=input_noisy[11,2],h3=input_noisy[11,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[12],c(Kd1=Kd1[12],Kd2=Kd2[12],Kd3=Kd3[12],k0=k0[12],ksyn=ksyn[12],kdeg=kdeg[12],h1=input_noisy[12,1],h2=input_noisy[12,2],h3=input_noisy[12,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[13],c(Kd1=Kd1[13],Kd2=Kd2[13],Kd3=Kd3[13],k0=k0[13],ksyn=ksyn[13],kdeg=kdeg[13],h1=input_noisy[13,1],h2=input_noisy[13,2],h3=input_noisy[13,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[14],c(Kd1=Kd1[14],Kd2=Kd2[14],Kd3=Kd3[14],k0=k0[14],ksyn=ksyn[14],kdeg=kdeg[14],h1=input_noisy[14,1],h2=input_noisy[14,2],h3=input_noisy[14,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[15],c(Kd1=Kd1[15],Kd2=Kd2[15],Kd3=Kd3[15],k0=k0[15],ksyn=ksyn[15],kdeg=kdeg[15],h1=input_noisy[15,1],h2=input_noisy[15,2],h3=input_noisy[15,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[16],c(Kd1=Kd1[16],Kd2=Kd2[16],Kd3=Kd3[16],k0=k0[16],ksyn=ksyn[16],kdeg=kdeg[16],h1=input_noisy[16,1],h2=input_noisy[16,2],h3=input_noisy[16,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[17],c(Kd1=Kd1[17],Kd2=Kd2[17],Kd3=Kd3[17],k0=k0[17],ksyn=ksyn[17],kdeg=kdeg[17],h1=input_noisy[17,1],h2=input_noisy[17,2],h3=input_noisy[17,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[18],c(Kd1=Kd1[18],Kd2=Kd2[18],Kd3=Kd3[18],k0=k0[18],ksyn=ksyn[18],kdeg=kdeg[18],h1=input_noisy[18,1],h2=input_noisy[18,2],h3=input_noisy[18,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[19],c(Kd1=Kd1[19],Kd2=Kd2[19],Kd3=Kd3[19],k0=k0[19],ksyn=ksyn[19],kdeg=kdeg[19],h1=input_noisy[19,1],h2=input_noisy[19,2],h3=input_noisy[19,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[20],c(Kd1=Kd1[20],Kd2=Kd2[20],Kd3=Kd3[20],k0=k0[20],ksyn=ksyn[20],kdeg=kdeg[20],h1=input_noisy[20,1],h2=input_noisy[20,2],h3=input_noisy[20,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[21],c(Kd1=Kd1[21],Kd2=Kd2[21],Kd3=Kd3[21],k0=k0[21],ksyn=ksyn[21],kdeg=kdeg[21],h1=input_noisy[21,1],h2=input_noisy[21,2],h3=input_noisy[21,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[22],c(Kd1=Kd1[22],Kd2=Kd2[22],Kd3=Kd3[22],k0=k0[22],ksyn=ksyn[22],kdeg=kdeg[22],h1=input_noisy[22,1],h2=input_noisy[22,2],h3=input_noisy[22,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[23],c(Kd1=Kd1[23],Kd2=Kd2[23],Kd3=Kd3[23],k0=k0[23],ksyn=ksyn[23],kdeg=kdeg[23],h1=input_noisy[23,1],h2=input_noisy[23,2],h3=input_noisy[23,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[24],c(Kd1=Kd1[24],Kd2=Kd2[24],Kd3=Kd3[24],k0=k0[24],ksyn=ksyn[24],kdeg=kdeg[24],h1=input_noisy[24,1],h2=input_noisy[24,2],h3=input_noisy[24,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[25],c(Kd1=Kd1[25],Kd2=Kd2[25],Kd3=Kd3[25],k0=k0[25],ksyn=ksyn[25],kdeg=kdeg[25],h1=input_noisy[25,1],h2=input_noisy[25,2],h3=input_noisy[25,3])),
                      get(paste("P1_GRN",GRN_ind,sep=""))(R_basal[26],c(Kd1=Kd1[26],Kd2=Kd2[26],Kd3=Kd3[26],k0=k0[26],ksyn=ksyn[26],kdeg=kdeg[26],h1=input_noisy[26,1],h2=input_noisy[26,2],h3=input_noisy[26,3])))
      GRS_all<-rbind(GRS_all,Gene_exp_rep)
    }
  }
  
  write(rep_ind,file="Nondyn_noise_generate_count.txt",append = T)
  return(GRS_all)
}
# save(dynamic_exp_rep,file="dynamic_exp_rep2.RData")
# save(dynamic_exp_rep,file="dynamic_exp_rep.RData")
# load("dynamic_exp_rep.RData")
logic_name<-c("1a2a3","1o2o3","1o(2a3)","1a(2o3)","2o(1a3)","2a(1o3)","3o(1a2)","3a(1o2)")
para_name<-paste(expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,1],
                 expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,2],
                 expand.grid(c("S","M","W"),c("S","M","W"),c("S","M","W"))[,3],sep="")[-27]
GRS_name<-paste(rep(logic_name,each=26),rep(para_name,8))
rownames(dynamic_exp_true)<-GRS_name
###############Save ground truth data
Nonredundent_name<-read.table(file = "Nonredundent_name.txt",sep = "\t")

Nondyn_active_full<-dynamic_exp_true[which(rownames(dynamic_exp_true)%in%apply(Nonredundent_name,1,as.character)),]
Nondyn_active_bio<-(100/Nondyn_active_full[,601])*Nondyn_active_full
Nondyn_WT_bio<-Nondyn_active_full[,601]
Nondyn_active_GT<-Nondyn_active_bio[,(c(0,15,30,60)/int_length+1)+rep(seq(0,801*25,801),each=4)]
Nondyn_active_GT_Full<-Nondyn_active_bio[,(c(0,0,2,14,15,16,29,30,31,58,60,60)/int_length+1)+rep(seq(0,801*25,801),each=12)]
#####Calclulate Slope from ground truth
ind_stdv<-rep(1:104,each=3)
slope_list_GT<-vector()
for(slope_delta in 1:60){
  ind_str<-c(0,0,0,(15-slope_delta),15,(15+slope_delta),(30-slope_delta),30,(30+slope_delta),(60-slope_delta),60,(60+slope_delta))
  ind_str[ind_str<0]<-0
  ind_str[ind_str>80]<-80
  
  tmp_length<-tapply(ind_str,rep(1:4,each=3),function(x) (x[3]-x[1]))
  tmp_length[1]<-1
  
  ind_temp<-(ind_str/int_length+1)+rep(seq(0,801*25,801),each=12)
  Nondyn_active_bio_f<-Nondyn_active_bio[,ind_temp]
  slope_GT<-vector()
  #####check
  for(i in 1:93){
    slope_GT<-c(slope_GT,tapply(Nondyn_active_bio_f[i,]/rep(rep(tmp_length,26),each=3),ind_stdv,function(x) (x[3]-x[1])))
  }
  slope_list_GT<-rbind(slope_list_GT,slope_GT)
}
########Save the data
save(Nondyn_WT_bio,file="Nondyn_WT_bio.RData")
save(Nondyn_active_bio,file="Nondyn_active_bio.RData")
save(Nondyn_active_GT,file="Nondyn_active_GT.RData")
save(slope_list_GT,file="slope_list_GT.RData")
save(Nondyn_active_GT_Full,file="Nondyn_active_GT_Full.RData")
# load("Nondyn_active_bio.RData")
# load("Nondyn_active_GT.RData")
# load("Nondyn_WT_bio.RData")
# load("slope_list_GT.RData")
########selected nonredundent activated noisy data
Nondyn_rep_M1.05<-list()
for(i in 1:1000){
  Nondyn_rep_M1.05[[i]]<-(100/Nondyn_WT_bio)*dynamic_exp_rep[[i]][which(rownames(dynamic_exp_true)%in%apply(Nonredundent_name,1,as.character)),]
  print(i)
}

save(Nondyn_rep_M1.05,file="Nondyn_rep_M1.05.RData")
###########################variance analysis
#########variance of input and parameter varied data
var_empirical_bio_para_M1.05<-matrix(0,ncol=104,nrow = 93)
ind_j<-rep(c(0,15,30,60),26)/int_length+rep(seq(0,801*25,801),each=4)+1

for(i in 1:93){
  for(j in 1:104){
    var_empirical_bio_para_M1.05[i,j]<-var(unlist(lapply(Nondyn_rep_M1.05,function(x) x[i,ind_j[j]])))
    # print(c(i,j))
  }
}

save(var_empirical_bio_para_M1.05,file="var_empirical_bio_para_M1.05.RData")
########################################
#Generate tempral variance: biological + technical
########################################
# Nondyn_active_bio<-Nondyn_rep_M1.05[[1]]
temp_noise_generate<-function(bio_var,tech_var,Nondyn_active_bio,rep_num,random_ind){
  int_length<-0.1
  Nondyn_active_temp<-vector()
  temp_ind_f_all<-vector()
  for(rep_ind in 1:rep_num){
    set.seed(rep_ind+random_ind)
    temp_bio_var<-rnorm(104*93, mean = 0, sd = bio_var)
    temp_tech_var<-rnorm(104*93, mean = 0, sd = tech_var)
    temp_var<-temp_tech_var+temp_bio_var
    temp_ind<-round(temp_var,digits=1)+rep(c(0,15,30,60),26*93)
    
    temp_ind[temp_ind<0]<-0
    temp_ind[temp_ind>80]<-80
    temp_ind[seq(1,9672,4)]<-0
    
    temp_ind_f<-t(apply(t(matrix(((temp_ind/int_length)+1),nrow=104)),1,function(x) x+rep(seq(0,801*25,801),each=4)))
    temp_ind_f_all<-rbind(temp_ind_f_all,temp_ind_f)
    Nondyn_active_temp_each<-vector()
    for(i in 1:93){
      Nondyn_active_temp_each<-rbind(Nondyn_active_temp_each,Nondyn_active_bio[i,temp_ind_f[i,]])
    }
    Nondyn_active_temp<-rbind(Nondyn_active_temp,Nondyn_active_temp_each)
  }
  return(Nondyn_active_temp)
}
###########genrate temporal
Nondyn_active_temp_BT50_parV_rep1000_N_M1.05<-vector()
for(i in 1:rep_size){
  Nondyn_active_temp_BT50_parV_rep1000_N_M1.05<-rbind(Nondyn_active_temp_BT50_parV_rep1000_N_M1.05,temp_noise_generate(1.414*5,1.414*5,Nondyn_rep_M1.05[[i]],1,i))
  print(i)
}

rownames(Nondyn_active_temp_BT50_parV_rep1000_N_M1.05)<-rep(apply(Nonredundent_name,1,as.character),rep_size)
save(Nondyn_active_temp_BT50_parV_rep1000_N_M1.05,file="Nondyn_active_temp_BT50_parV_rep1000_N_M1.05.RData")
#########variance of simulated data
var_ij<-matrix(0,ncol=104,nrow = 93)
for(i in 1:93){
  for(j in 1:104){
    var_ij[i,j]<-var(Nondyn_active_temp_BT50_parV_rep1000_N_M1.05[(c(i)+seq(0,(rep_size-1)*93,93)),j])
  }
}
var_empirical_pure_bio_V0.01_BT50_rep1000_M1.05<-var_ij
save(var_empirical_pure_bio_V0.01_BT50_rep1000_M1.05,file="var_empirical_pure_bio_V0.01_BT50_rep1000_M1.05.RData")
#########Plot empirical-theoretical variance 
# load("var_empirical_pure_bio_V0.01_BT50_rep1000_M1.05.RData")
slope_ind<-round(2*2*deltaN*(70)^0.5)
if(slope_ind>60) slope_ind<-60
if(slope_ind<1) slope_ind<-1

pdf("Temp_var_estimate_BT50_parV_rep1000_M1.05.pdf",width=5,height=5)
plot(as.matrix(t(var_empirical_pure_bio_V0.01_BT50_rep1000_M1.05),ncol=1),as.matrix(slope_list_GT[slope_ind,]^2*70*deltaN^2,ncol=1),xlim=c(0,2000),ylim=c(0,2000),xlab="Variance of simulated data",ylab="Estimated variance from theory",pch=16)
abline(coef = c(0,1),col="red",lwd=2)
graphics.off()
########################################
#Generate basal level:0.05
########################################Uniform[0,1]
Nondyn_active_temp_BT50_parV_rep1000_N_M1.05[Nondyn_active_temp_BT50_parV_rep1000_N_M1.05==0]<-0.05
########################################
#Generate sequencing variance: var~mean^2
########################################
##############select the data
############################################Noise generate
set.seed(1)
########################noise generation function
# dyn_data<-Nondyn_active_temp
# delta<-0.02
# var_g_t<-40
# sample_num<-1000
noise_generate<-function(dyn_data,delta,sample_num,random_ind){
  set.seed(random_ind)
  dyn_data_Noise<-vector()
  Var_list<-vector()
  for(GRS_index in 1:dim(dyn_data)[1]){
    var_RNA<-delta*dyn_data[GRS_index,]^2
    dyn_data_list<-vector()
    for(i in 1:104){
      if(dyn_data[GRS_index,i]>3){
        data_sample<-rnorm(sample_num*10,mean=dyn_data[GRS_index,i],sd=sqrt(var_RNA[i]))
        while(any(abs(subset(data_sample,data_sample>0)[1:sample_num]-dyn_data[GRS_index,i])>5*sqrt(var_RNA[i]))){
          data_sample<-rnorm(sample_num*10,mean=dyn_data[GRS_index,i],sd=sqrt(var_RNA[i]))
        }
        dyn_data_list<-cbind(dyn_data_list,subset(data_sample,data_sample>0)[1:sample_num])
      } else{
        data_sample<-rgamma(sample_num*10,shape=(dyn_data[GRS_index,i])^2/var_RNA[i],scale = var_RNA[i]/(dyn_data[GRS_index,i]))
        while(any(abs(subset(data_sample,data_sample>0)[1:sample_num]-dyn_data[GRS_index,i])>5*sqrt(var_RNA[i]))){
          data_sample<-rnorm(sample_num*10,mean=dyn_data[GRS_index,i],sd=sqrt(var_RNA[i]))
        }
        dyn_data_list<-cbind(dyn_data_list,subset(data_sample,data_sample>0)[1:sample_num])
      }
    }
    dyn_data_Noise<-rbind(dyn_data_Noise,as.vector(t(dyn_data_list)))
    Var_list<-rbind(Var_list,var_RNA)
  }
  rownames(dyn_data_Noise)<-rownames(Nondyn_active_bio)
  return(list(dyn_data_Noise=dyn_data_Noise,Var_list=Var_list))
}
################1000################
##replicates: TB5,TT5
rep_num<-1000
Nondyn_bio_V0.01_BT50_rep1000_M1.05<-vector()
for(i in 1:rep_num){
  Nondyn_bio_V0.01_BT50_rep1000_M1.05<-cbind(Nondyn_bio_V0.01_BT50_rep1000_M1.05,noise_generate(Nondyn_active_temp_BT50_parV_rep1000_N_M1.05[((93*i-92):(93*i)),],0.01*deltaN^2,1,i)$dyn_data_Noise)
  print(i)
}
##########
save(Nondyn_bio_V0.01_BT50_rep1000_M1.05,file="Nondyn_bio_V0.01_BT50_rep1000_M1.05.RData")
################################################################
################################################################
################Empirical Variance################
#############TB5_TT5
var_empirical_bio_V0.01_BT50_rep1000_M1.05<-matrix(0,ncol=104,nrow = 93)
for(i in 1:93){
  for(j in 1:104){
    var_empirical_bio_V0.01_BT50_rep1000_M1.05[i,j]<-var(Nondyn_bio_V0.01_BT50_rep1000_M1.05[i,(j+seq(0,104*999,104))])
  }
}
save(var_empirical_bio_V0.01_BT50_rep1000_M1.05,file="var_empirical_bio_V0.01_BT50_rep1000_M1.05.RData")
########################################
################################Plot dot
########################################
# pdf("Dyn_Noise_bio_V0.01_BT50_rep1000.pdf",width=8,height=8)
# par(mfrow=c(3,3),mai = c(0.3, 0.3, 0.3, 0.3))
# for(plot_input in 1:9){
#   plot(c(0,15,30,60),Nondyn_active_GT[80,(plot_input*4-3):(plot_input*4)],col="black",xlab='',ylab='',type="p",panel.first = grid(),cex.lab=1.3,pch=16,ylim=c(0,200),main=paste("TF Perturb",plot_input))
#   for(i in 1:1000){
#     points(c(0,15,30,60),Nondyn_bio_V0.01_BT50_rep1000[80,((plot_input*4-3):(plot_input*4)+104*(i-1))],col="grey",pch=16)
#   }
#   points(c(0,15,30,60),Nondyn_active_GT[80,(plot_input*4-3):(plot_input*4)],col="black",pch=16)
# }
# graphics.off()

########################################
################################Plot input
########################################
noise_ind<-1
input_noisy_list_all<-vector()
input_per_list_all<-vector()
para_noisy_list_all<-vector()
para_per_list_all<-vector()

for(rep_ind in 1:2){
  input_noisy_list<-vector()
  para_noisy_list<-vector()
  input_per_list<-vector()
  para_per_list<-vector()
  
  for(GRN_ind in 1:8){
    for(para_ind in 1:dim(para_list)[1]){
      input_noisy_list<-rbind(input_noisy_list,as.vector(t(input_raw*(1+input_error[,(noise_ind*3-2):(noise_ind*3)]))))
      
      temp<-(1+input_error[,(noise_ind*3-2):(noise_ind*3)])
      temp[input_raw==0]<-0
      
      input_per_list<-rbind(input_per_list,as.vector(t(temp)))
      
      Kd1<-para_list[para_ind,1]*exp(para_noise[(26*(5*noise_ind-4)-25):(26*(5*noise_ind-4))])
      Kd2<-para_list[para_ind,2]*exp(para_noise[(26*(5*noise_ind-3)-25):(26*(5*noise_ind-3))])
      Kd3<-para_list[para_ind,3]*exp(para_noise[(26*(5*noise_ind-2)-25):(26*(5*noise_ind-2))])
      ksyn<-1*exp(para_noise[(26*(5*noise_ind-1)-25):(26*(5*noise_ind-1))])
      kdeg<-0.1*exp(para_noise[(130*noise_ind-25):(130*noise_ind)])
      k0<-0+k0_noise[(26*noise_ind-25):(26*noise_ind)]
      
      para_noisy_list<-rbind(para_noisy_list,as.vector(t(cbind(Kd1,Kd2,Kd3,ksyn,kdeg,k0))))
      para_per_list<-rbind(para_per_list,as.vector(t(cbind(exp(para_noise[(26*(5*noise_ind-4)-25):(26*(5*noise_ind-4))]),
                                                           exp(para_noise[(26*(5*noise_ind-3)-25):(26*(5*noise_ind-3))]),
                                                           exp(para_noise[(26*(5*noise_ind-2)-25):(26*(5*noise_ind-2))]),
                                                           exp(para_noise[(26*(5*noise_ind-1)-25):(26*(5*noise_ind-1))]),
                                                           exp(para_noise[(130*noise_ind-25):(130*noise_ind)]),
                                                           k0_noise[(26*noise_ind-25):(26*noise_ind)]))))
      
      
      noise_ind<-noise_ind+1
    }
  }
  if(length(input_noisy_list_all)==0) {
    input_noisy_list_all<-input_noisy_list
    para_noisy_list_all<-para_noisy_list
    
    input_per_list_all<-input_per_list
    para_per_list_all<-para_per_list
  }  else{
    input_noisy_list_all<-cbind(input_noisy_list_all,input_noisy_list)
    para_noisy_list_all<-cbind(para_noisy_list_all,para_noisy_list)
    
    input_per_list_all<-cbind(input_per_list_all,input_per_list)
    para_per_list_all<-cbind(para_per_list_all,para_per_list)
  }
}
rownames(input_noisy_list_all)<-GRS_name
rownames(para_noisy_list_all)<-GRS_name
rownames(input_per_list_all)<-GRS_name
rownames(para_per_list_all)<-GRS_name

input_plot_raw1<-input_noisy_list_all[which(rownames(input_noisy_list_all)%in%apply(Nonredundent_name,1,as.character)),]
input_plot_raw2<-rbind(input_plot_raw1[,1:78],input_plot_raw1[,79:156])
input_plot<-input_plot_raw2[rep(c(1:93),each=2)+rep(c(0,93),93),]

para_plot_raw1<-para_noisy_list_all[which(rownames(para_noisy_list_all)%in%apply(Nonredundent_name,1,as.character)),]
para_plot_raw2<-rbind(para_plot_raw1[,1:156],para_plot_raw1[,157:312])
para_plot<-para_plot_raw2[rep(c(1:93),each=2)+rep(c(0,93),93),]
para_plot[,seq(4,156,6)]<-rep((100/Nondyn_WT_bio),each=2)*para_plot[,seq(4,156,6)]

input_per_plot_raw1<-input_per_list_all[which(rownames(input_per_list_all)%in%apply(Nonredundent_name,1,as.character)),]
input_per_plot_raw2<-rbind(input_per_plot_raw1[,1:78],input_per_plot_raw1[,79:156])
input_per_plot<-input_per_plot_raw2[rep(c(1:93),each=2)+rep(c(0,93),93),]

para_per_plot_raw1<-para_per_list_all[which(rownames(para_per_list_all)%in%apply(Nonredundent_name,1,as.character)),]
para_per_plot_raw2<-rbind(para_per_plot_raw1[,1:156],para_per_plot_raw1[,157:312])
para_per_plot<-para_per_plot_raw2[rep(c(1:93),each=2)+rep(c(0,93),93),]
para_per_plot[,seq(6,156,6)]<-NA
#################################Plot##############
pdf("Input_noise_plot_M1.05.pdf",width=10,height=14)
pheatmap(input_plot,cellheight=3, cellwidth = 5,gaps_col = seq(3,75,3),gaps_row = seq(2,93*2,2),legend=T,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()
##############################################
para_plot[para_plot<0.02]<-0.02
pdf("Para_noise_plot_M1.05.pdf",width=20,height=14)
pheatmap(log2(para_plot),cellheight=3, cellwidth = 5,gaps_col = seq(6,156,6),gaps_row = seq(2,93*2,2),legend=T,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()
##############################################
pdf("Input_per_plot_M1.05.pdf",width=10,height=14)
pheatmap(input_per_plot,cellheight=3, cellwidth = 5,gaps_col = seq(3,75,3),gaps_row = seq(2,93*2,2),legend=T,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()
##############################################
pdf("Para_per_plot_M1.05.pdf",width=20,height=14)
pheatmap(para_per_plot,cellheight=3, cellwidth = 5,gaps_col = seq(6,156,6),gaps_row = seq(2,93*2,2),legend=T,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()
##############################################
# load("Nondyn_bio_V0.01_BT50_rep1000_M1.05.RData")
gene_raw<-rbind(Nondyn_bio_V0.01_BT50_rep1000_M1.05[,1:104],Nondyn_bio_V0.01_BT50_rep1000_M1.05[,105:208])
gene_plot<-gene_raw[rep(c(1:93),each=2)+rep(c(0,93),93),]

pdf("Gene_plot_M1.05.pdf",width=20,height=14)
pheatmap(gene_plot,cellheight=3, cellwidth = 5,gaps_col = seq(4,104,4),gaps_row = seq(2,93*2,2),legend=T,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()

gene_raw<-rbind(Nondyn_bio_V0.01_BT50_rep1000_M1.05[,1:104],Nondyn_bio_V0.01_BT50_rep1000_M1.05[,105:208])
Nondyn_active_GT_N<-gene_raw/Nondyn_active_GT[rep(1:93,2),]
Nondyn_active_GT_N[Nondyn_active_GT_N==Inf]<-NA

pdf("Gene_per_M1.05.pdf",width=20,height=14)
pheatmap(Nondyn_active_GT_N,cellheight=3, cellwidth = 5,gaps_col = seq(4,104,4),gaps_row = seq(2,93*2,2),legend=T,
         show_rownames = F,show_colnames = F,cluster_cols = F,cluster_rows = F)
graphics.off()
######################################
################################################################################
load("var_empirical_bio_V0.01_BT50_rep1000_M1.05.RData")
Average_dis<-summary(as.vector(var_empirical_bio_V0.01_BT50_rep1000_M1.05))[4]
pdf("Average_dis_V0.01_BT50.pdf",width=6,height=6)
hist(as.vector(var_empirical_bio_V0.01_BT50_rep1000_M1.05),breaks=30,xlab = c("Variance"),main=paste("Mean of the distribution",round(Average_dis,digit=2)),cex=1.5)
graphics.off()
######################################
source("Error_Model_Infer_BT50_M1.05.R")
######################################
# source("/home/ning/Error_Model/Optimizaton_upgrade/Nondyn_opt_BT50_0.02/Nondyn_opt_BT50_M1.05_up.R")
# source("/home/ning/Error_Model/Optimizaton_upgrade/Nondyn_Con_opt_BT50_M1.05/Nondyn_Con_opt_BT50_M1.05_up.R")