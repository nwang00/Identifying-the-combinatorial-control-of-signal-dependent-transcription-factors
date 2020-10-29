library(compiler)
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
###################################################
P1_GRN1_fit_raw<-function(t, state, parameters) {
  
  t_sudo<-t*(1/int_length)
  t_f<-floor(t_sudo)+1
  t_a_f<-t_sudo-(t_f-1)
  
  G <- (1-k0)*((t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,])*(t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,])*(t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,])/((parameters["Kd1"]+t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,])*(parameters["Kd2"]+t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,])*(parameters["Kd3"]+t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,])))+k0
  # ptint(c(G,state))
  
  dRNA <-parameters["ksyn"]*G-parameters["kdeg"]*state
  
  list(dRNA)
}


P1_GRN2_fit_raw<-function(t, state, parameters) {
  
  t_sudo<-t*(1/int_length)
  t_f<-floor(t_sudo)+1
  t_a_f<-t_sudo-(t_f-1)
  
  G<-(1-k0)*((1-parameters["Kd1"]*parameters["Kd2"]*parameters["Kd3"]/((parameters["Kd1"]+(t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,]))*(parameters["Kd2"]+(t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,]))*(parameters["Kd3"]+(t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,])))))+k0
  
  dRNA <-parameters["ksyn"]*G-parameters["kdeg"]*state
  
  list(dRNA)
}

P1_GRN3_fit_raw<-function(t, state, parameters) {
  
  t_sudo<-t*(1/int_length)
  t_f<-floor(t_sudo)+1
  t_a_f<-t_sudo-(t_f-1)
  
  G<-(1-k0)*(1-(parameters["Kd1"]/(parameters["Kd1"]+t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,]))*(1-(t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,])*(t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,])/((parameters["Kd2"]+t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,])*(parameters["Kd3"]+t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,]))))+k0
  
  dRNA <-parameters["ksyn"]*G-parameters["kdeg"]*state
  
  list(dRNA)
}

P1_GRN4_fit_raw<-function(t, state, parameters) {
  
  t_sudo<-t*(1/int_length)
  t_f<-floor(t_sudo)+1
  t_a_f<-t_sudo-(t_f-1)
  
  G<-(1-k0)*(((t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,])/(parameters["Kd1"]+t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,]))*(1-parameters["Kd2"]*parameters["Kd3"]/((parameters["Kd2"]+(t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,]))*(parameters["Kd3"]+(t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,])))))+k0
  
  dRNA <-parameters["ksyn"]*G-parameters["kdeg"]*state
  
  list(dRNA)
}

P1_GRN5_fit_raw<-function(t, state, parameters) {
  
  t_sudo<-t*(1/int_length)
  t_f<-floor(t_sudo)+1
  t_a_f<-t_sudo-(t_f-1)
  
  G<-(1-k0)*(1-(parameters["Kd2"]/(parameters["Kd2"]+t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,]))*(1-(t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,])*(t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,])/((parameters["Kd1"]+t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,])*(parameters["Kd3"]+t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,]))))+k0
  
  dRNA <-parameters["ksyn"]*G-parameters["kdeg"]*state
  
  list(dRNA)
}

P1_GRN6_fit_raw<-function(t, state, parameters) {
  
  t_sudo<-t*(1/int_length)
  t_f<-floor(t_sudo)+1
  t_a_f<-t_sudo-(t_f-1)
  
  G<-(1-k0)*(((t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,])/(parameters["Kd2"]+t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,]))*(1-parameters["Kd1"]*parameters["Kd3"]/((parameters["Kd1"]+(t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,]))*(parameters["Kd3"]+(t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,])))))+k0
  
  dRNA <-parameters["ksyn"]*G-parameters["kdeg"]*state
  
  list(dRNA)
}

P1_GRN7_fit_raw<-function(t, state, parameters) {
  
  t_sudo<-t*(1/int_length)
  t_f<-floor(t_sudo)+1
  t_a_f<-t_sudo-(t_f-1)
  
  G<-(1-k0)*(1-(parameters["Kd3"]/(parameters["Kd3"]+t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,]))*(1-(t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,])*(t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,])/((parameters["Kd1"]+t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,])*(parameters["Kd2"]+t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,]))))+k0
  
  dRNA <-parameters["ksyn"]*G-parameters["kdeg"]*state
  
  list(dRNA)
}

P1_GRN8_fit_raw<-function(t, state, parameters) {
  
  t_sudo<-t*(1/int_length)
  t_f<-floor(t_sudo)+1
  t_a_f<-t_sudo-(t_f-1)
  
  G<-(1-k0)*(((t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,])/(parameters["Kd3"]+t_a_f*cJun[t_f+1,]+(1-t_a_f)*cJun[t_f,]))*(1-parameters["Kd1"]*parameters["Kd2"]/((parameters["Kd1"]+(t_a_f*NFkB[t_f+1,]+(1-t_a_f)*NFkB[t_f,]))*(parameters["Kd2"]+(t_a_f*IRF3[t_f+1,]+(1-t_a_f)*IRF3[t_f,])))))+k0
  
  dRNA <-parameters["ksyn"]*G-parameters["kdeg"]*state
  
  list(dRNA)
}

P1_GRN1_fit <- cmpfun(P1_GRN1_fit_raw)
P1_GRN2_fit <- cmpfun(P1_GRN2_fit_raw)
P1_GRN3_fit <- cmpfun(P1_GRN3_fit_raw)
P1_GRN4_fit <- cmpfun(P1_GRN4_fit_raw)
P1_GRN5_fit <- cmpfun(P1_GRN5_fit_raw)
P1_GRN6_fit <- cmpfun(P1_GRN6_fit_raw)
P1_GRN7_fit <- cmpfun(P1_GRN7_fit_raw)
P1_GRN8_fit <- cmpfun(P1_GRN8_fit_raw)
##########################################################
##########################################################
##########################################################
P1_fun_GRNs<-function(X){
  Kd1<-10^(X[1])
  Kd2<-10^(X[2])
  Kd3<-10^(X[3])
  ksyn<-10^(X[4])
  kdeg<-10^(X[5])
  k0<-k0
  

  R_basal<-ksyn*(get(paste("G_basal",GRN_ind,sep=""))(Kd1,Kd2,Kd3,k0,Basal_level))/kdeg

  Total_mse<-target_fun(matrix(ode(y=c(RNA1=R_basal,RNA2=R_basal,RNA3=R_basal,RNA4=R_basal,RNA5=R_basal), 
                                   times = caRNA_time,func =get(paste("P1_GRN",GRN_ind,"_fit",sep="")),parms = c(Kd1=Kd1,Kd2=Kd2,Kd3=Kd3,ksyn=ksyn,kdeg=kdeg,k0=k0),hmax=1)[,-1],ncol=1))

  if(is.infinite(Total_mse)){return(10^8)}
  if(is.na(Total_mse)){return(10^8)}
  return(Total_mse)
}
