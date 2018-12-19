#name = "conbase"
#name = "monovar"
#name = "sccaller"

#flat = T


plot_simulation_stats_each_ds <- function(name,flat=T){

if (flat){
   infile = sprintf("../data/flat_dp/stats_all4_sim_%s.csv",name)
   outfile = sprintf("../plots/%s_flat_dp.pdf",name)
   savefile = sprintf("../data/flat_dp/stat_matrix_%s.Rdata",name)    
}else {
   infile = sprintf("../data/sim_dp/stats_all4_sim_%s.csv",name)
   outfile = sprintf("../plots/%s_sim_dp.pdf",name)
   savefile = sprintf("../data/sim_dp/stat_matrix_%s.Rdata",name)    
}

D = read.table(infile, header=T, sep=",", row.names = 1)

pdf(outfile, width=14,height=10)
###########################################
# define colors:

# T-het      red	
# T-hom      darkgreen	
# F-het      darkmagenta	
# F-hom      darkblue  	
# NA-het     darkgray - grey10
# NA-hom     darkgray - grey30   
# T-het_wEAL lightred	  
# T-hom_wEAL lightgreen	  
# F-het_wEAL ligthmagenta	  
# F-hom_wEAL libhtblue  	  
# NA-het_wEAL    lightgray - grey 50
# NA-hom_wEAL    lightgray - grey30    


library(RColorBrewer)
spec <- brewer.pal(11,"Spectral")
make_lighter <- function(color, factor=0.5){
  x = rgb2hsv(col2rgb(color))
  v = pmax(pmin(x[3, ] + factor, 1), 0)
  hsv(h = x[1, ], s = x[2, ], v = v)
  }


coldef <- c(spec[8],spec[10],spec[2], spec[4], "grey30","grey40")
coldef <- c(coldef, sapply(coldef, make_lighter, factor=1))
coldef[11:12] <- c("grey50","grey60")

#coldef <- c("red","darkgreen","magenta4","darkblue","grey20","grey30","pink","lightgreen","magenta1","lightblue","grey50","grey70")
names(coldef) <- rownames(D)[1:12]



pl = 1:12
###############################################
# plot for all with fSNV = 1 and eal=0
eal0 = which(D["fEAL",]==0)
ado0 = which(D["fADO",]==0)



par(mfrow=c(1,2), xpd=NA, oma=c(1,0,0,5), cex=0.7)
barplot(as.matrix(D[pl,eal0]), names.arg = D["fADO",eal0], main=sprintf("%s fEAL=0",name), col=coldef, las=2)
mtext("fADO",1)
barplot(as.matrix(D[pl,ado0]), names.arg = D["fEAL",ado0], main=sprintf("%s fADO=0",name), col=coldef, las=2)
legend(12,6000,fill=coldef, legend=names(coldef), bty='n', cex=1)
mtext("fEAL",1)



###############################################
# calculate TP etc. 
# merge data from w_EAL and without



make_stats <- function(M){
# true false calls:
TC <- 1:2
FC <- 3:4

M2 <- rbind(M[1:6,]+M[7:12,], M[21:24,])

# number of GT calls.
M2["nGT",] = (colSums(M2[TC,])+colSums(M2[FC,]))

# fraction of correct genotype tGT among all GT calls
M2["fTrueGT",] = colSums(M2[TC,])/M2["nGT",]

# number of deteced M2ut calls.
M2["nMut",] <- colSums(M2[c("T-het","F-hom"),])

# fraction of correct mut calls among mut calls
M2["fTrueMut",] = colSums(M2["T-het",]) / M2["nMut",]

# total simulated mutations per sample - should be nMut + NA-het
M2["simMut",] <- M2["nMut",] + M2["NA-het",]

# fraction of correct mut calls among all true muts
M2["fDMut",] = colSums(M2["T-het",]) / M2["simMut",]

return(M2)
}

# summary:
# - nGT - total number of genotype calls.
# - fTrueGT - fraction of all GT-calls that has the correct genotype (both 1 and 0)
# - nMut - total calls that predict a mutation 
# - fTrueMut - fraction of Mut calls (has mutation) that are correct)
# - simMut - total simulated mutations across all loci in all cells.
# - fDMut - fraction of true mutations that are detected.

D2 <- make_stats(D)

###############################################
# Make matrices with ADO and EAL on X/Y axis 
# number of detected sites
# number of detected mutations
# fraction of correct GT-calls vs total calls.

library(ggplot2)
library(gridExtra)


# transpose the dfs
D2 <- data.frame(t(D2))



nogrid <- theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
scale <- scale_colour_gradientn(colours = c("yellow","red","black"), limits = c(0,1), na.value = "gray")
bx <- scale_x_continuous(breaks = seq(0,0.9,by=0.2), limits = c(0,0.9))
by <- scale_y_continuous(breaks = seq(0,0.9,by=0.2), limits = c(0,0.9))


p1 <- ggplot(D2, aes(x=fADO, y=fEAL, colour = fTrueMut)) + geom_point(data=D2[D2$nMut > 0,],aes(size=nMut)) + ggtitle(sprintf("%s - Found Mut / Pred Mut",name)) + scale + nogrid + bx + by
p2 <- ggplot(D2, aes(x=fADO, y=fEAL, colour = fTrueGT)) + geom_point(data=D2[D2$nGT > 0,],aes(size=nGT)) + ggtitle(sprintf("%s - Correct GT / nGT",name)) + scale + nogrid + bx + by
p3 <- ggplot(D2, aes(x=fADO, y=fEAL, colour = fDMut)) + geom_point(data=D2[D2$simMut > 0,],aes(size=simMut)) + ggtitle(sprintf("%s - Found Mut / Simulated Muts",name)) + scale + nogrid + bx + by


grid.arrange(p1,p2,p3, ncol=2)


# summary:
# - nGT - total number of genotype calls.
# - fTrueGT - fraction of all GT-calls that has the correct genotype (both 1 and 0)
# - nMut - total calls that predict a mutation 
# - fTrueMut - fraction of Mut calls (has mutation) that are correct
# - simMut - total simulated mutations across all loci in all cells.
# - fDMut - fraction of simulated mutations that are detected.

###############################################
# FP / FN

# TP -> Mut_2
# FP -> Mut_2_w_EAL
# Simulated positives: nSim -> nSite-nEAL
# Simulated negatives: nEAL 

# TN -> nEAL-Mut_2_wEAL
# FN -> nSim-Mut_2

# Sensitivity -> TP/ nSim (total simulated positives) 
# Specificity -> TN / nEAL (total simulated negatives)

# TPR = TP / (TP + FN)
# FPR = FP / (FP + TN) == 1-Spec.
# 


nSite <- 306
#coldef.tp <- c("green","red","magenta","gray")
coldef.tp <- c(spec[8],spec[10],spec[2], "grey40")

names(coldef.tp) <- c("TP","TN","FP","FN")

# for cutoff 2
D["nSim",] <- nSite - D["nEAL",]
D["TP2",] <- D["Mut_2",]
D["TN2",] <- D["nEAL",]-D["Mut_2_wEAL",]
D["FP2",] <- D["Mut_2_wEAL",]
D["FN2",] <- D["nSim",] - D["Mut_2",]
D["Sens2",] <- D["TP2",] / D["nSim",]
D["Spec2",] <- D["TN2",] / D["nEAL",]
D["Det2",] <- D["FP2",] + D["TP2",]
D["FDR2",] <- D["FP2",] / D["Det2",]

# for cutoff 5
D["TP5",] <- D["Mut_5",]
D["TN5",] <- D["nEAL",]-D["Mut_5_wEAL",]
D["FP5",] <- D["Mut_5_wEAL",]
D["FN5",] <- D["nSim",] - D["Mut_5",]
D["Sens5",] <- D["TP5",] / D["nSim",]
D["Spec5",] <- D["TN5",] / D["nEAL",]

# for cutoff 10
D["TP10",] <- D["Mut_10",]
D["TN10",] <- D["nEAL",]-D["Mut_10_wEAL",]
D["FP10",] <- D["Mut_10_wEAL",]
D["FN10",] <- D["nSim",] - D["Mut_10",]
D["Sens10",] <- D["TP10",] / D["nSim",]
D["Spec10",] <- D["TN10",] / D["nEAL",]



# with eal=0 the Spec becomes NA - 
D[is.na(D)] <- -1


# make barplot with eal=0/ado0
# cutoff 2
pl2 <- 26:29

par(mfrow=c(2,2), xpd=NA, oma=c(1,0,1,5), cex=0.7)
barplot(as.matrix(D[pl2,eal0]), names.arg = D["fADO",eal0], main=sprintf("%s fEAL=0",name), col=coldef.tp, las=2)
mtext("fADO",1)
barplot(as.matrix(D[pl2,ado0]), names.arg = D["fEAL",ado0], main=sprintf("%s fADO=0",name), col=coldef.tp, las=2)
legend(12,300,fill=coldef.tp, legend=names(coldef.tp), bty='n', cex=1.2, title="cutoff 2")
mtext("fEAL",1)


# cutoff 5
pl2 <- 34:38
barplot(as.matrix(D[pl2,eal0]), names.arg = D["fADO",eal0], main=sprintf("%s fEAL=0",name), col=coldef.tp, las=2)
mtext("fADO",1)
barplot(as.matrix(D[pl2,ado0]), names.arg = D["fEAL",ado0], main=sprintf("%s fADO=0",name), col=coldef.tp, las=2)
legend(12,300,fill=coldef.tp, legend=names(coldef.tp), bty='n', cex=1.2, title="cutoff 5")
mtext("fEAL",1)


# cutoff 10
pl2 <- 40:43
barplot(as.matrix(D[pl2,eal0]), names.arg = D["fADO",eal0], main=sprintf("%s fEAL=0",name), col=coldef.tp, las=2)
mtext("fADO",1)
barplot(as.matrix(D[pl2,ado0]), names.arg = D["fEAL",ado0], main=sprintf("%s fADO=0",name), col=coldef.tp, las=2)
legend(12,300,fill=coldef.tp, legend=names(coldef.tp), bty='n', cex=1.2, title="cutoff 10")
mtext("fEAL",1)


# make sens & spec plot
par(mfrow = c(2,2))
plot(as.matrix(D["fADO",eal0]),as.matrix(D["Sens2",eal0]), type='b', ylab="Sensitivity",xlab="fADO")
plot(as.matrix(D["fADO",eal0]),as.matrix(D["Spec2",eal0]), type='b', ylab="Specificity",xlab="fADO")

plot(as.matrix(D["fEAL",ado0]),as.matrix(D["Sens2",ado0]), type='b', ylab="Sensitivity",xlab="fEAL")
plot(as.matrix(D["fEAL",ado0]),as.matrix(D["Spec2",ado0]), type='b', ylab="Specificity",xlab="fEAL")
mtext("Cutoff2",3,outer=T)


values = seq(0,0.9,by=0.1)
cols = rainbow(10)



par(mfrow=c(2,2))
# cutoff2
plot(as.matrix(D["fADO",eal0]),as.matrix(D["Sens2",eal0]), type='l', ylab="Sensitivity",xlab="fADO", col=cols[1],ylim=c(0,1))
for (i in 2:10){
sel = which(D["fEAL",]==values[i])
    lines(as.matrix(D["fADO",sel]),as.matrix(D["Sens2",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fEAL", cex=0.7, bty='n')
mtext("cutoff 2",3)

# cutoff5
plot(as.matrix(D["fADO",eal0]),as.matrix(D["Sens5",eal0]), type='l', ylab="Sensitivity",xlab="fADO", col=cols[1],ylim=c(0,1))
for (i in 2:10){
sel = which(D["fEAL",]==values[i])
    lines(as.matrix(D["fADO",sel]),as.matrix(D["Sens5",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fEAL", cex=0.7, bty='n')
mtext("cutoff 5",3)

# cutoff10
plot(as.matrix(D["fADO",eal0]),as.matrix(D["Sens10",eal0]), type='l', ylab="Sensitivity",xlab="fADO", col=cols[1],ylim=c(0,1))
for (i in 2:10){
sel = which(D["fEAL",]==values[i])
    lines(as.matrix(D["fADO",sel]),as.matrix(D["Sens10",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fEAL",cex=0.7, bty='n')
mtext("cutoff 10",3)



# same for ado0

par(mfrow=c(2,2))
# cutoff2
plot(as.matrix(D["fEAL",ado0]),as.matrix(D["Sens2",ado0]), type='l', ylab="Sensitivity",xlab="fEAL", col=cols[1],ylim=c(0,1))
for (i in 2:10){
sel = which(D["fADO",]==values[i])
    lines(as.matrix(D["fEAL",sel]),as.matrix(D["Sens2",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fADO", cex=0.7, bty='n')
mtext("cutoff 2",3)


# cutoff5
plot(as.matrix(D["fEAL",ado0]),as.matrix(D["Sens5",ado0]), type='l', ylab="Sensitivity",xlab="fEAL", col=cols[1], ylim=c(0,1))
for (i in 2:10){
sel = which(D["fADO",]==values[i])
    lines(as.matrix(D["fEAL",sel]),as.matrix(D["Sens5",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fADO", cex=0.7, bty='n')
mtext("cutoff 5",3)

# cutoff10
plot(as.matrix(D["fEAL",ado0]),as.matrix(D["Sens10",ado0]), type='l', ylab="Sensitivity",xlab="fEAL", col=cols[1], ylim=c(0,1))
for (i in 2:10){
sel = which(D["fADO",]==values[i])
    lines(as.matrix(D["fEAL",sel]),as.matrix(D["Sens10",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fADO", cex=0.7, bty='n')
mtext("cutoff 10",3)


# same type of plot for specificity
par(mfrow=c(2,2))
# cutoff2
plot(as.matrix(D["fADO",eal0]),as.matrix(D["Spec2",eal0]), type='l', ylab="Specificity",xlab="fADO", col=cols[1],ylim=c(0,1))
for (i in 2:10){
sel = which(D["fEAL",]==values[i])
    lines(as.matrix(D["fADO",sel]),as.matrix(D["Spec2",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fADO", cex=0.7, bty='n')
mtext("cutoff 2",3)


# cutoff5
plot(as.matrix(D["fADO",eal0]),as.matrix(D["Spec5",eal0]), type='l', ylab="Specificity",xlab="fADO", col=cols[1],ylim=c(0,1))
for (i in 2:10){
sel = which(D["fEAL",]==values[i])
    lines(as.matrix(D["fADO",sel]),as.matrix(D["Spec5",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fADO", cex=0.7, bty='n')
mtext("cutoff 5",3)


# cutoff10
plot(as.matrix(D["fADO",eal0]),as.matrix(D["Spec10",eal0]), type='l', ylab="Specificity",xlab="fADO", col=cols[1],ylim=c(0,1))
for (i in 2:10){
sel = which(D["fEAL",]==values[i])
    lines(as.matrix(D["fADO",sel]),as.matrix(D["Spec10",sel]), col=cols[i])
}
legend("bottomleft",as.character(values), fill=cols, title="fADO", cex=0.7, bty='n')
mtext("cutoff 10",3)




# also do dotplots for sens spec, FPR
D2 <- data.frame(t(D))


p1 <- ggplot(D2, aes(x=fADO, y=fEAL, colour = Sens2)) + geom_point(data=D2[D2$nSim > 0,],aes(size=nSim)) + ggtitle(sprintf("%s - Sensitivity cut2",name)) + scale + nogrid + bx + by
p2 <- ggplot(D2, aes(x=fADO, y=fEAL, colour = Spec2)) + geom_point(data=D2[D2$nSim > 0,],aes(size=nSim)) + ggtitle(sprintf("%s - Specificity cut2",name)) + scale + nogrid + bx + by

p3 <- ggplot(D2, aes(x=fADO, y=fEAL, colour = FDR2)) + geom_point(data=D2[D2$Det2 > 0,],aes(size=Det2)) + ggtitle(sprintf("%s - False Discovery Rate cut2",name)) + scale + nogrid + bx + by

grid.arrange(p1,p2,p3,ncol=2)

# ROC plots
par(mfrow=c(2,2))
plot(1-as.matrix(D["Spec2",]),as.matrix(D["Sens2",]), ylab="Sensitivity", xlab = "1-Specificity", main="All settings")
eal0.5 <- which(D["fEAL",] == 0.5)
plot(1-as.matrix(D["Spec2",eal0.5]),as.matrix(D["Sens2",eal0.5]), ylab="Sensitivity", xlab = "1-Specificity", main="eal=0.5")
ado0.5  <- which(D["fADO",] == 0.5)
plot(1-as.matrix(D["Spec2",ado0.5]),as.matrix(D["Sens2",ado0.5]), ylab="Sensitivity", xlab = "1-Specificity", main="ado=0.5")

dev.off()
save(D2,file=savefile)

}


# ROC - FPR (x) and TPR (y)
# TPR == Sensitivity
# FPR = 1-Spec
# TNR = Spec = TN/N = TN/(TN+FP)

