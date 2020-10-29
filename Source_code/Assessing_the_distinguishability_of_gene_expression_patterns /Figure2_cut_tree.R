##############################################################
##############################################################
rm(list=ls())
library("ggplot2")
library(dendextend)
library(RColorBrewer)
load("Identifiability_7sets_SSE_basalN.RData")
load("Identifiability_all_SSE_basalN.RData")
load("Identifiability_SSE_basalN_Delayed.RData")
load("Identifiability_SSE_basalN_Gradient.RData")
load("Identifiability_SSE_basalN_Transient.RData")

con_reorder<-c(25,21,9,19,7,3,1)
con_ind_f<-unlist(lapply(con_reorder,function(x) seq((x*4-3),x*4,1)))
gene_exp<-dynamic_norm[,con_ind_f]
colnames(SSE_7sets)<-rownames(gene_exp)
rownames(SSE_7sets)<-rownames(gene_exp)
SSE_7sets[SSE_7sets<0.1]<-0

hclust7_single <-hclust(as.dist(SSE_7sets), method = "single", members = NULL)
dend7 <- hclust7_single %>% as.dendrogram
##########################################################################
##########################################################################
##########################################################################
SSE_all[SSE_all<0.1]<-0

hclustall_single<-hclust(as.dist(SSE_all), method = "single", members = NULL)
dendall <-hclustall_single%>%as.dendrogram
hclust_27HML<-hclustall_single
##########################################################################
SSE_all_delayed[SSE_all_delayed<0.1]<-0

hclustdelayed_single<-hclust(as.dist(SSE_all_delayed), method = "single", members = NULL)
dend_delayed<-hclustdelayed_single%>%as.dendrogram
##########################################################################
SSE_all_gradient[SSE_all_gradient<0.1]<-0

hclustgradient_single<-hclust(as.dist(SSE_all_gradient), method = "single", members = NULL)
dend_gradient<-hclustgradient_single%>%as.dendrogram
##########################################################################
SSE_all_transient[SSE_all_transient<0.1]<-0

hclusttransient_single<-hclust(as.dist(SSE_all_transient), method = "single", members = NULL)
dend_transient<-hclusttransient_single%>%as.dendrogram
##########################################################################
##########################################################################threshold-cluster number plot
##########################################################################
duplicated_id_7<-which(duplicated(names(heights_per_k.dendrogram(dend7))))
name_list_7<-heights_per_k.dendrogram(dend7)
if(length(duplicated_id_7)!=0){
  name_list_nonR_7<-name_list_7[-duplicated_id_7]
}
name_list_nonR_7<-name_list_7
name_list_nonR_7<-name_list_nonR_7[-which(name_list_nonR_7<0)]

den_plot_7<-name_list_nonR_7
den_plot_select_log10_7<-round(log10(den_plot_7),digits = 2)[-1]

plot_thre_7<-cbind(seq(1,5,0.01),rep(0,length(seq(1,5,0.01))))

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
name_list_nonR_all<-name_list_all[-duplicated_id_all]

den_plot_all<-name_list_nonR_all
den_plot_select_log10_all<-round(log10(den_plot_all),digits = 2)[-1]

plot_thre_all<-cbind(seq(1,5,0.01),rep(0,length(seq(1,5,0.01))))

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
if(length(duplicated_id_delayed)!=0){
  name_list_nonR_delayed<-name_list_delayed[-duplicated_id_delayed]
}
name_list_nonR_delayed<-name_list_delayed
name_list_nonR_delayed<-name_list_nonR_delayed[-which(name_list_nonR_delayed<0)]

den_plot_delayed<-name_list_nonR_delayed
den_plot_select_log10_delayed<-round(log10(den_plot_delayed),digits = 2)[-1]

plot_thre_delayed<-cbind(seq(1,5,0.01),rep(0,length(seq(1,5,0.01))))

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
  name_list_nonR_gradient<-name_list_nonR_gradient[-which(name_list_nonR_gradient<0)]
}

den_plot_gradient<-name_list_nonR_gradient
den_plot_select_log10_gradient<-round(log10(den_plot_gradient),digits = 2)[-1]

plot_thre_gradient<-cbind(seq(1,5,0.01),rep(0,length(seq(1,5,0.01))))

for(i in 1:length(den_plot_select_log10_gradient)){
  ind_s<-intersect(which(plot_thre_gradient[,1]<=den_plot_select_log10_gradient[i]),which(plot_thre_gradient[,1]>den_plot_select_log10_gradient[i+1]))
  plot_thre_gradient[ind_s,2]<-as.numeric(names(den_plot_select_log10_gradient)[i])
}
plot_thre_gradient[which(plot_thre_gradient[,1]>=den_plot_select_log10_gradient[1]),2]<-as.numeric(names(den_plot_select_log10_gradient)[1])
plot_thre_gradient[which(plot_thre_gradient[,1]<=tail(den_plot_select_log10_gradient,n=1)),2]<-as.numeric(tail(names(den_plot_select_log10_gradient),n=1))
##########################################################################
##########################################################################All transient
duplicated_id_transient<-which(duplicated(names(heights_per_k.dendrogram(dend_transient))))
name_list_transient<-heights_per_k.dendrogram(dend_transient)
if(length(duplicated_id_transient)!=0){
  name_list_nonR_transient<-name_list_transient[-duplicated_id_transient]
}
name_list_nonR_transient<-name_list_transient
name_list_nonR_transient<-name_list_nonR_transient[-which(name_list_nonR_transient<0)]

den_plot_transient<-name_list_nonR_transient
den_plot_select_log10_transient<-round(log10(den_plot_transient),digits = 2)[-1]

plot_thre_transient<-cbind(seq(1,5,0.01),rep(0,length(seq(1,5,0.01))))

for(i in 1:length(den_plot_select_log10_transient)){
  ind_s<-intersect(which(plot_thre_transient[,1]<=den_plot_select_log10_transient[i]),which(plot_thre_transient[,1]>den_plot_select_log10_transient[i+1]))
  plot_thre_transient[ind_s,2]<-as.numeric(names(den_plot_select_log10_transient)[i])
}
plot_thre_transient[which(plot_thre_transient[,1]>=den_plot_select_log10_transient[1]),2]<-as.numeric(names(den_plot_select_log10_transient)[1])
plot_thre_transient[which(plot_thre_transient[,1]<=tail(den_plot_select_log10_transient,n=1)),2]<-as.numeric(tail(names(den_plot_select_log10_transient),n=1))
##########################
pdf("Cluster_treeheight_single_MDTG_revised.pdf",width=6,height=6)
plot(plot_thre_7[,1],plot_thre_7[,2],ylim=c(0,110),col="black",xlab="",ylab="",type="o",cex=1.5,cex.axis=1.5,pch=16,lwd=3)
lines(plot_thre_all[,1],plot_thre_all[,2],col=brewer.pal(4, "Set1")[1],type="o",pch=16,cex=1.5,lwd=3)
lines(plot_thre_delayed[,1],plot_thre_delayed[,2],col=brewer.pal(4, "Set1")[2],type="o",pch=16,cex=1.5,lwd=3)
lines(plot_thre_gradient[,1],plot_thre_gradient[,2],col=brewer.pal(4, "Set1")[3],type="o",pch=16,cex=1.5,lwd=3)
lines(plot_thre_transient[,1],plot_thre_transient[,2],col=brewer.pal(4, "Set1")[4],type="o",pch=16,cex=1.5,lwd=3)
# legend("topright",legend=c("7 Datasets","Amplitute","Delayed","Gradient","Transient"),
#        col=c("black",brewer.pal(4, "Set1")),pch=c(16,16,16,16,16), cex=1)
graphics.off()