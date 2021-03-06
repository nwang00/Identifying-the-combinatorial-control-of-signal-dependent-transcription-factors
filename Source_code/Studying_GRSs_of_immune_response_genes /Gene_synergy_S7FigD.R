load("/home/diane/Clicc_Share/Ning_paper/Likelihood_data.RData")

LL_list <- cbind(LL_list, data.frame("1" = NA, "2" = NA, "3" = NA, "1o2" = NA, "1o3" = NA, "2o3" = NA, "1a2" = NA, "1a3" = NA, "2a3" = NA, "0" = NA, check.names = F))

names(LL_list) <- gsub(" ", "", names(LL_list))

for (gate in names(LL_list)[1:8]){
  strength_gate <- RS_matrix[grep(pattern = gate, row.names(RS_matrix), fixed=T),]
  genes <- sapply(row.names(strength_gate), function(x){strsplit(x,' ')[[1]][1]})
  
  print(all(genes == row.names(LL_list)))
  
  LL <- LL_list[,gate]
  if( gate == "1a2a3" ){
    simplified_genes <- apply(strength_gate,1,function(x){any(x == "N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "0"] <- apply(cbind(LL_list[simplified_genes, "0"],LL[simplified_genes]), 1, min, na.rm=T)
  }else if( gate == "1o2o3" ){
    simplified_genes <- apply(strength_gate,1,function(x){all(x == "N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "0"] <-  apply(cbind(LL_list[simplified_genes, "0"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){all(x[c(1,2)] == "N") & x[3]!="N"})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "3"] <- apply(cbind(LL_list[simplified_genes, "3"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){all(x[c(1,3)] == "N") & x[2]!="N"})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "2"] <- apply(cbind(LL_list[simplified_genes, "2"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){all(x[c(2,3)] == "N") & x[1]!="N"})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1"] <- apply(cbind(LL_list[simplified_genes, "1"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[1] == "N" & all(x[c(2,3)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "2o3"] <- apply(cbind(LL_list[simplified_genes, "2o3"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[2] == "N" & all(x[c(1,3)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1o3"] <- apply(cbind(LL_list[simplified_genes, "1o3"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[3] == "N" & all(x[c(1,2)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1o2"] <- apply(cbind(LL_list[simplified_genes, "1o2"],LL[simplified_genes]), 1, min, na.rm=T)
  }else if( gate == "1o(2a3)" ){
    simplified_genes <- apply(strength_gate,1,function(x){x[1] == "N" & any(x[c(2,3)] == "N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "0"] <-  apply(cbind(LL_list[simplified_genes, "0"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){any(x[c(2,3)] == "N") & x[1]!="N"})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1"] <-  apply(cbind(LL_list[simplified_genes, "1"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[1] == "N" & all(x[c(2,3)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "2a3"] <- apply(cbind(LL_list[simplified_genes, "2a3"],LL[simplified_genes]), 1, min, na.rm=T)
  }else if( gate == "1a(2o3)" ){
    simplified_genes <- apply(strength_gate,1,function(x){x[1] == "N" | all(x[c(2,3)] == "N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "0"] <-  apply(cbind(LL_list[simplified_genes, "0"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[2] == "N" & all(x[c(1,3)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1a3"] <- apply(cbind(LL_list[simplified_genes, "1a3"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[3] == "N" & all(x[c(1,2)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1a2"] <- apply(cbind(LL_list[simplified_genes, "1a2"],LL[simplified_genes]), 1, min, na.rm=T)
  }else if( gate == "2o(1a3)" ){
    simplified_genes <- apply(strength_gate,1,function(x){x[2] == "N" & any(x[c(1,3)] == "N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "0"] <-  apply(cbind(LL_list[simplified_genes, "0"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){any(x[c(1,3)] == "N") & x[2]!="N"})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "2"] <-  apply(cbind(LL_list[simplified_genes, "2"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[2] == "N" & all(x[c(1,3)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1a3"] <- apply(cbind(LL_list[simplified_genes, "1a3"],LL[simplified_genes]), 1, min, na.rm=T)
  }else if( gate == "2a(1o3)" ){
    simplified_genes <- apply(strength_gate,1,function(x){x[2] == "N" | all(x[c(1,3)] == "N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "0"] <-  apply(cbind(LL_list[simplified_genes, "0"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[1] == "N" & all(x[c(2,3)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "2a3"] <- apply(cbind(LL_list[simplified_genes, "2a3"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[3] == "N" & all(x[c(1,2)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1a2"] <- apply(cbind(LL_list[simplified_genes, "1a2"],LL[simplified_genes]), 1, min, na.rm=T)
  }else if( gate == "3o(1a2)" ){
    simplified_genes <- apply(strength_gate,1,function(x){x[3] == "N" & any(x[c(1,2)] == "N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "0"] <-  apply(cbind(LL_list[simplified_genes, "0"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){any(x[c(1,2)] == "N") & x[3]!="N"})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "3"] <-  apply(cbind(LL_list[simplified_genes, "3"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[3] == "N" & all(x[c(1,2)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1a2"] <- apply(cbind(LL_list[simplified_genes, "1a2"],LL[simplified_genes]), 1, min, na.rm=T)
  }else if( gate == "3a(1o2)" ){
    simplified_genes <- apply(strength_gate,1,function(x){x[3] == "N" || all(x[c(1,2)] == "N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "0"] <-  apply(cbind(LL_list[simplified_genes, "0"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[1] == "N" & all(x[c(2,3)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "2a3"] <- apply(cbind(LL_list[simplified_genes, "2a3"],LL[simplified_genes]), 1, min, na.rm=T)
    
    simplified_genes <- apply(strength_gate,1,function(x){x[2] == "N" & all(x[c(1,3)]!="N")})
    LL_list[simplified_genes, gate] <- NA
    LL_list[simplified_genes, "1a3"] <- apply(cbind(LL_list[simplified_genes, "1a3"],LL[simplified_genes]), 1, min, na.rm=T)
  }
}

# # plot 1 - 2
# or <- apply(LL_list[,c(2,3,5, 8, 12)], 1, min, na.rm=T)
# and <- apply(LL_list[,c(1,4,6,7, 15)], 1, min, na.rm=T)
# best <-  apply(LL_list, 1, function(x){idx <- which.min(x); idx %in% c(1:8,12,15)})
# 
# dat <- data.frame(or = or, and = and, best=best)
# p <- ggplot(dat, aes(x=or,y=and, col=best)) + geom_point()
# p <- p + geom_abline(slope=1,)
# p
# 
# 
# or <- apply(LL_list[,c(2,3,5,8,12,9:11,13:14,16:17)], 1, min, na.rm=T)
# and <- apply(LL_list[,c(1,4,6,7,15)], 1, min, na.rm=T)
# 
# dat <- data.frame(or = or, and = and)
# p <- ggplot(dat, aes(x=or,y=and, col=best)) + geom_point()
# p <- p + geom_abline(slope=1,)
# p
# 
# # plot 1 - 3 
# or <- LL_list[,c(13)]
# and <- LL_list[,c(16)]
# best <-  apply(LL_list, 1, function(x){idx <- which.min(x); idx %in% c(1:8,13,16)})
# 
# dat <- data.frame(or = or, and = and, best=best)
# p <- ggplot(dat, aes(x=or,y=and, col=best)) + geom_point()
# p <- p + geom_abline(slope=1)
# p
# 
# or <- apply(LL_list[,c(2,3,6,7,13)], 1, min, na.rm=T)
# and <- apply(LL_list[,c(1,4,5,8,16)], 1, min, na.rm=T)
# best <-  apply(LL_list, 1, function(x){idx <- which.min(x); idx %in% c(1:8,13,16)})
# 
# dat <- data.frame(or = or, and = and, best=best)
# p <- ggplot(dat, aes(x=or,y=and, col=best)) + geom_point()
# p <- p + geom_abline(slope=1,)
# p
# 
# or <- apply(LL_list[,c(2,3,6,7,13,9:12,14:15,17)], 1, min, na.rm=T)
# and <- apply(LL_list[,c(1,4,5,8,16)], 1, min, na.rm=T)
# best <-  apply(LL_list, 1, function(x){idx <- which.min(x); idx %in% c(1:8,13,16)})
# 
# dat <- data.frame(or = or, and = and)
# p <- ggplot(dat, aes(x=or,y=and)) + geom_point()
# p <- p + geom_abline(slope=1,)
# p
# 
# # plot 2 - 3 
# or <- apply(LL_list[,c(2,4,5,7,14)], 1, min, na.rm=T)
# and <- apply(LL_list[,c(1,3,6,8,17)], 1, min, na.rm=T)
# best <-  apply(LL_list, 1, function(x){idx <- which.min(x); idx %in% c(1:8,14,17)})
# 
# dat <- data.frame(or = or, and = and, best=best)
# p <- ggplot(dat, aes(x=or,y=and, col=best)) + geom_point()
# p <- p + geom_abline(slope=1,)
# p
# 
# or <- apply(LL_list[,c(2,4,5,7,14,9:13,15:16)], 1, min, na.rm=T)
# and <- apply(LL_list[,c(1,3,6,8,17)], 1, min, na.rm=T)
# 
# dat <- data.frame(or = or, and = and, best=best)
# p <- ggplot(dat, aes(x=or,y=and)) + geom_point()
# p <- p + geom_abline(slope=1,)
# p
# 
# 
# # or vs and pure logic
# or <- apply(LL_list[,c(2,9:14)], 1, min, na.rm=T)
# and <- apply(LL_list[,c(1,15:17)], 1, min, na.rm=T)
# 
# dat <- data.frame(or = or, and = and)
# p <- ggplot(dat, aes(x=or,y=and))
# p <- p + geom_polygon(data = data.frame(x=c(-20, 200, -20), y=c(-20, 200, 200)), aes(x,y),fill="cyan", alpha=0.2)
# p <- p + geom_polygon(data = data.frame(x=c(-20, 200, 200), y=c(-20, 200, -20)), aes(x,y),fill="orange", alpha=0.2)
# p <- p + geom_point()
# p <- p + geom_abline(slope=1, col="red",linetype=2)
# p <- p + labs(x="No synergy", y = "Synergy")
# p <- p + theme_bw() + theme(aspect.ratio = 1)
# p <- p + annotate(geom="text", x=0, y=80, hjust=0, label="More likely no synergy")
# p <- p + annotate(geom="text", x=100, y=0, hjust=1, label="More likely synergy")
# p <- p + coord_cartesian(xlim = c(-10,120), ylim=c(-10,120))
# p

# any synergy vs no synergy
or <- apply(LL_list[,c(2,9:14,18)], 1, min, na.rm=T)
and <- apply(LL_list[,c(1,3:8,15:17)], 1, min, na.rm=T)

dat <- data.frame(or = or, and = and)
p <- ggplot(dat, aes(x=or,y=and))
p <- p + geom_polygon(data = data.frame(x=c(-20, 200, -20), y=c(-20, 200, 200)), aes(x,y),fill="cyan", alpha=0.2)
p <- p + geom_polygon(data = data.frame(x=c(-20, 200, 200), y=c(-20, 200, -20)), aes(x,y),fill="orange", alpha=0.2)
p <- p + geom_point()
p <- p + geom_abline(slope=1, col="red",linetype=2)
p <- p + labs(x="No synergy", y = "Synergy")
p <- p + theme_bw() + theme(aspect.ratio = 1)
p <- p + annotate(geom="text", x=0, y=80, hjust=0, label="More likely no synergy")
p <- p + annotate(geom="text", x=100, y=0, hjust=1, label="More likely synergy")
p <- p + coord_cartesian(xlim = c(-10,120), ylim=c(-10,120))
p

