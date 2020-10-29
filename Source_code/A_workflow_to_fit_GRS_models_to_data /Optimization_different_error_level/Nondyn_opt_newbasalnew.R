####################
P1_GRN1<-function(R_basal, parameters) {
  G<-(1-parameters["k0"])*(parameters["h1"]*parameters["h2"]*parameters["h3"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd2"]+parameters["h2"])*(parameters["Kd3"]+parameters["h3"])))+parameters["k0"]
  return(G*parameters["ksyn"]/parameters["kdeg"]+(R_basal-G*parameters["ksyn"]/parameters["kdeg"])*exp(-parameters["kdeg"]*caRNA_time))
}
P1_GRN2<-function(R_basal, parameters) {
  G<-(1-parameters["k0"])*(1-parameters["Kd1"]*parameters["Kd2"]*parameters["Kd3"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd2"]+parameters["h2"])*(parameters["Kd3"]+parameters["h3"])))+parameters["k0"]
  return(G*parameters["ksyn"]/parameters["kdeg"]+(R_basal-G*parameters["ksyn"]/parameters["kdeg"])*exp(-parameters["kdeg"]*caRNA_time))
}
P1_GRN3<-function(R_basal, parameters) {
  G<-(1-parameters["k0"])*(1-(parameters["Kd1"]/(parameters["Kd1"]+parameters["h1"]))*(1-parameters["h2"]*parameters["h3"]/((parameters["Kd2"]+parameters["h2"])*(parameters["Kd3"]+parameters["h3"]))))+parameters["k0"]
  return(G*parameters["ksyn"]/parameters["kdeg"]+(R_basal-G*parameters["ksyn"]/parameters["kdeg"])*exp(-parameters["kdeg"]*caRNA_time))
}
P1_GRN4<-function(R_basal, parameters) {
  G<-(1-parameters["k0"])*((parameters["h1"]/(parameters["Kd1"]+parameters["h1"]))*(1-parameters["Kd2"]*parameters["Kd3"]/((parameters["Kd2"]+parameters["h2"])*(parameters["Kd3"]+parameters["h3"]))))+parameters["k0"]
  return(G*parameters["ksyn"]/parameters["kdeg"]+(R_basal-G*parameters["ksyn"]/parameters["kdeg"])*exp(-parameters["kdeg"]*caRNA_time))
}
P1_GRN5<-function(R_basal, parameters) {
  G<-(1-parameters["k0"])*(1-(parameters["Kd2"]/(parameters["Kd2"]+parameters["h2"]))*(1-parameters["h1"]*parameters["h3"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd3"]+parameters["h3"]))))+parameters["k0"]
  return(G*parameters["ksyn"]/parameters["kdeg"]+(R_basal-G*parameters["ksyn"]/parameters["kdeg"])*exp(-parameters["kdeg"]*caRNA_time))
}
P1_GRN6<-function(R_basal, parameters) {
  G<-(1-parameters["k0"])*((parameters["h2"]/(parameters["Kd2"]+parameters["h2"]))*(1-parameters["Kd1"]*parameters["Kd3"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd3"]+parameters["h3"]))))+parameters["k0"]
  return(G*parameters["ksyn"]/parameters["kdeg"]+(R_basal-G*parameters["ksyn"]/parameters["kdeg"])*exp(-parameters["kdeg"]*caRNA_time))
}
P1_GRN7<-function(R_basal, parameters) {
  G<-(1-parameters["k0"])*(1-(parameters["Kd3"]/(parameters["Kd3"]+parameters["h3"]))*(1-parameters["h1"]*parameters["h2"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd2"]+parameters["h2"]))))+parameters["k0"]
  return(G*parameters["ksyn"]/parameters["kdeg"]+(R_basal-G*parameters["ksyn"]/parameters["kdeg"])*exp(-parameters["kdeg"]*caRNA_time))
}
P1_GRN8<-function(R_basal, parameters) {
  G<-(1-parameters["k0"])*((parameters["h3"]/(parameters["Kd3"]+parameters["h3"]))*(1-parameters["Kd1"]*parameters["Kd2"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd2"]+parameters["h2"]))))+parameters["k0"]
  return(G*parameters["ksyn"]/parameters["kdeg"]+(R_basal-G*parameters["ksyn"]/parameters["kdeg"])*exp(-parameters["kdeg"]*caRNA_time))
}
################################
##################################
##################################
G_GRN1<-function(parameters) {
  G<-(1-parameters["k0"])*(parameters["h1"]*parameters["h2"]*parameters["h3"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd2"]+parameters["h2"])*(parameters["Kd3"]+parameters["h3"])))+parameters["k0"]
  return(G)
}

G_GRN2<-function(parameters) {
  G<-(1-parameters["k0"])*(1-parameters["Kd1"]*parameters["Kd2"]*parameters["Kd3"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd2"]+parameters["h2"])*(parameters["Kd3"]+parameters["h3"])))+parameters["k0"]
  return(G)
}

G_GRN3<-function(parameters) {
  G<-(1-parameters["k0"])*(1-(parameters["Kd1"]/(parameters["Kd1"]+parameters["h1"]))*(1-parameters["h2"]*parameters["h3"]/((parameters["Kd2"]+parameters["h2"])*(parameters["Kd3"]+parameters["h3"]))))+parameters["k0"]
  return(G)
}

G_GRN4<-function( parameters) {
  G<-(1-parameters["k0"])*((parameters["h1"]/(parameters["Kd1"]+parameters["h1"]))*(1-parameters["Kd2"]*parameters["Kd3"]/((parameters["Kd2"]+parameters["h2"])*(parameters["Kd3"]+parameters["h3"]))))+parameters["k0"]
  return(G)
}

G_GRN5<-function( parameters) {
  G<-(1-parameters["k0"])*(1-(parameters["Kd2"]/(parameters["Kd2"]+parameters["h2"]))*(1-parameters["h1"]*parameters["h3"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd3"]+parameters["h3"]))))+parameters["k0"]
  return(G)
}

G_GRN6<-function( parameters) {
  G<-(1-parameters["k0"])*((parameters["h2"]/(parameters["Kd2"]+parameters["h2"]))*(1-parameters["Kd1"]*parameters["Kd3"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd3"]+parameters["h3"]))))+parameters["k0"]
  return(G)
}

G_GRN7<-function( parameters) {
  G<-(1-parameters["k0"])*(1-(parameters["Kd3"]/(parameters["Kd3"]+parameters["h3"]))*(1-parameters["h1"]*parameters["h2"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd2"]+parameters["h2"]))))+parameters["k0"]
  return(G)
}

G_GRN8<-function( parameters) {
  G<-(1-parameters["k0"])*((parameters["h3"]/(parameters["Kd3"]+parameters["h3"]))*(1-parameters["Kd1"]*parameters["Kd2"]/((parameters["Kd1"]+parameters["h1"])*(parameters["Kd2"]+parameters["h2"]))))+parameters["k0"]
  return(G)
}
##################################
##################################
##################################
##################################
##################################
##################################
G_basal1<-function(Kd1,Kd2,Kd3,k0,Basal_level){
  return((1-k0)*((Basal_level)^3/((Kd1+Basal_level)*(Kd2+Basal_level)*(Kd3+Basal_level)))+k0)
}

G_basal2<-function(Kd1,Kd2,Kd3,k0,Basal_level){
  return((1-k0)*((1-Kd1*Kd2*Kd3/((Kd1+Basal_level)*(Kd2+Basal_level)*(Kd3+Basal_level))))+k0)
}

G_basal3<-function(Kd1,Kd2,Kd3,k0,Basal_level){
  return((1-k0)*(1-(Kd1/(Kd1+Basal_level))*(1-Basal_level^2/((Kd2+Basal_level)*(Kd3+Basal_level))))+k0)
}

G_basal4<-function(Kd1,Kd2,Kd3,k0,Basal_level){
  return((1-k0)*((Basal_level/(Kd1+Basal_level))*(1-Kd2*Kd3/((Kd2+Basal_level)*(Kd3+Basal_level))))+k0)
}

G_basal5<-function(Kd1,Kd2,Kd3,k0,Basal_level){
  return((1-k0)*(1-(Kd2/(Kd2+Basal_level))*(1-Basal_level^2/((Kd1+Basal_level)*(Kd3+Basal_level))))+k0)
}

G_basal6<-function(Kd1,Kd2,Kd3,k0,Basal_level){
  return((1-k0)*((Basal_level/(Kd2+Basal_level))*(1-Kd1*Kd3/((Kd1+Basal_level)*(Kd3+Basal_level))))+k0)
}

G_basal7<-function(Kd1,Kd2,Kd3,k0,Basal_level){
  return((1-k0)*(1-(Kd3/(Kd3+Basal_level))*(1-Basal_level^2/((Kd1+Basal_level)*(Kd2+Basal_level))))+k0)
}

G_basal8<-function(Kd1,Kd2,Kd3,k0,Basal_level){
  return((1-k0)*((Basal_level/(Kd3+Basal_level))*(1-Kd1*Kd2/((Kd1+Basal_level)*(Kd2+Basal_level))))+k0)
}

P1_fun_GRNs<-function(X){
  # print(X)
  if(any(X[1:3]>2)|any(X[1:3]<(-2))|(X[6]>(1))|(X[6]<(0))) {return(10^8)}
  Kd1<-10^(X[1])
  Kd2<-10^(X[2])
  Kd3<-10^(X[3])
  ksyn<-10^(X[4])
  kdeg<-10^(X[5])
  k0<-X[6]
  
  # Kd1<-10^(0.356590974)
  # Kd2<-10^(-1.352434174)
  # Kd3<-10^(-0.028680619)
  # ksyn<-10^(1.264826761)
  # kdeg<-10^(-1.022643040)
  # k0<-0.002599714 
  
  #correct
  # Kd1<-10^(-1)
  # Kd2<-10^(0)
  # Kd3<-10^(1)
  # ksyn<-10.55
  # kdeg<-10^(-1)
  # k0<-0
  # # #best_fit
  # Kd1<-10^(-1.62)
  # Kd2<-10^(-0.095)
  # Kd3<-10^(1.05)
  # ksyn<-10^(0.91)
  # kdeg<-10^(-1.08)
  # k0<-0.0003

  #40
  #correct
  # Kd1<-10^(-1)
  # Kd2<-10^(-1)
  # Kd3<-10^(0)
  # ksyn<-10.54852
  # kdeg<-10^(-1)
  # k0<-0
  # # #best_fit
  # Kd1<-10^(-1.3211072421)
  # Kd2<-10^(0.0555631110)
  # Kd3<-10^(-0.6961266731)
  # ksyn<-10^(0.9333445404)
  # kdeg<-10^(-1.0864844054)
  # k0<-0.0001743699
  # # 681.7661
  # 
  # #70
  # #correct
  # Kd1<-10^(-1)
  # Kd2<-10^(-1)
  # Kd3<-10^(0)
  # ksyn<-11.56018
  # kdeg<-10^(-1)
  # k0<-0
  # # #best_fit
  # Kd1<-10^(-0.3654665228)
  # Kd2<-10^(-0.7637120744)
  # Kd3<-10^(0.2200285179)
  # ksyn<-10^(0.9272315376)
  # kdeg<-10^(-1.3843749226)
  # k0<-0.0008822873
  # 549.4992
  
  
  parS<-c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,k0=k0,ksyn=ksyn,kdeg=kdeg)
  R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg
  # print(R_basal)
  if(is.nan(R_basal)|(R_basal<0)) {return(10^8)}
  # print("yes")
  Total_mse<-target_fun(c(get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[1,1],h2=input_list[1,2],h3=input_list[1,3])),
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
                          get(paste("P1_GRN",GRN_ind,sep=""))(R_basal,c(parS,h1=input_list[26,1],h2=input_list[26,2],h3=input_list[26,3]))))
  # print(Total_mse)
  if(is.infinite(Total_mse)){return(10^8)}
  if(is.na(Total_mse)){return(10^8)}
  return(Total_mse)
}
