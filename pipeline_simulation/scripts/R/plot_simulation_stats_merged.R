
flat = T
all <- list()

library(RColorBrewer)
spec <- brewer.pal(11,"Spectral")

names = c("conbase","monovar","sccaller","lira")
coldef2 = spec[c(1,9,4,11)]


for (name in names){
if (flat){
   savefile = sprintf("../data/flat_dp/stat_matrix_%s.Rdata",name)
   outfile <- "../plots/merged_plots_flat_dp.pdf"
}else {
   savefile = sprintf("../data/sim_dp/stat_matrix_%s.Rdata",name)
   outfile <- "../plots/merged_plots.pdf"
}
load(savefile)
all[[name]] <- D2
}


names(coldef2) <- names(all)


make_lighter <- function(color, factor=0.5){
  x = rgb2hsv(col2rgb(color))
  v = pmax(pmin(x[3, ] + factor, 1), 0)
  hsv(h = x[1, ], s = x[2, ], v = v)
  }


coldef <- c(spec[8],spec[10],spec[2], spec[4], "grey30","grey60")
names(coldef) <- colnames(D2)[1:6]


coldef.tp <- c(spec[8],spec[10],spec[2],"grey40")
names(coldef.tp) <- c("TP","TN","FP","FN")




pdf(outfile,width=14,height=12, useDingbats=F)

###############################################
# plot GT detection eal=0
###############################################
# plot for all with fSNV = 1 and eal=0
D <- t(all[[1]])
eal0 <- which(D["fEAL",]==0)
ado0 <- which(D["fADO",]==0)

pl = 1:6

par(mfrow=c(2,2), xpd=NA, oma=c(1,0,0,5), cex=0.7)
for (name in names){
    D <- t(all[[name]])
    barplot(as.matrix(D[pl,eal0]), names.arg = D["fADO",eal0], main=sprintf("%s fEAL=0",name), col=coldef, las=2)
    mtext("fADO",1)
}
legend(12,6000,fill=coldef, legend=names(coldef), bty='n', cex=1)



###############################################
# plot GT detection eal=0, but with all bars beside eachoter
###############################################

D.merge <- Reduce(rbind,lapply(all, function(x) x[eal0,pl]))
D.merge$Method <- rep(names(all),each=10)

if (length(all) == 3){
# reorder
o <- as.vector(sapply(seq(1,10), function(x) x+c(0,10,20)))
D.merge <- D.merge[o,]

sp <- rep(c(1.5,0.5,0.5),10)
#cc <- coldef2[D.merge[,7]]
n <- as.vector(sapply( D["fADO",eal0], function(x) paste(c("C","M","S"),x,sep="_")))
par(mfrow = c(1,1), xpd=NA, oma=c(1,0,0,5))
barplot(t(as.matrix(D.merge[,1:6])), col=coldef, space = sp, names.arg = n, las=2)
legend(55,6000,fill=coldef, legend=names(coldef), bty='n', cex=1)

}else {

o <- as.vector(sapply(seq(1,10), function(x) x+c(0,10,20,30)))
D.merge <- D.merge[o,]

sp <- rep(c(1.5,0.5,0.5,0.5),10)
#cc <- coldef2[D.merge[,7]]
n <- as.vector(sapply( D["fADO",eal0], function(x) paste(c("C","M","S","L"),x,sep="_")))
par(mfrow = c(1,1), xpd=NA, oma=c(1,0,0,5))
barplot(t(as.matrix(D.merge[,1:6])), col=coldef, space = sp, names.arg = n, las=2)
legend(70,6000,fill=coldef, legend=names(coldef), bty='n', cex=1)



}

###############################################
# plot TP/FP with eal0/ado0, cutoff = 2
###############################################

pl2 <- 26:29


par(mfrow=c(2,2), xpd=NA, oma=c(1,0,0,5), cex=0.7)
for (name in names){
    D <- t(all[[name]])
    barplot(as.matrix(D[pl2,eal0]), names.arg = D["fADO",eal0], main=sprintf("%s fEAL=0",name), col=coldef.tp, las=2)
    mtext("fADO",1)
    barplot(as.matrix(D[pl2,ado0]), names.arg = D["fEAL",ado0], main=sprintf("%s fADO=0",name), col=coldef.tp, las=2)
    mtext("fEAL",1)
    legend(12,300,fill=coldef.tp, legend=names(coldef.tp), bty='n', cex=1.2)
}


###############################################
# plot FDR vs Detection dotplots.
###############################################

library(ggplot2)
library(gridExtra)


nogrid <- theme_minimal() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
scale <- scale_colour_gradientn(colours = c("white","black"), limits = c(0,1), na.value = "gray")
#scale <- scale_colour_gradientn(colours = c("yellow","red","black"), limits = c(0,1), na.value = "gray")
bx <- scale_x_continuous(breaks = seq(0,0.9,by=0.2), limits = c(0,0.9))
by <- scale_y_continuous(breaks = seq(0,0.9,by=0.2), limits = c(0,0.9))

scale_size <- scale_size_identity(trans="sqrt",guide="legend")

plots <- list()
for (name in names){
    D2 <- all[[name]]
    plots[[name]] <- ggplot(D2, aes(x=fADO, y=fEAL,colour = FDR2)) + geom_point(data=D2[D2$Det2 > 0,],aes(size=Det2)) + geom_point(data=D2[D2$Det2 > 0,],aes(size=Det2),shape=1,colour="black") + ggtitle(sprintf("%s - False Discovery Rate",name)) + scale + nogrid + bx + by +  scale_size(limits = c(0,306),range=c(0,15))



}

grid.arrange(grobs = plots, ncol=2)


# Or do size by FP
plots <- list()
for (name in names){
    D2 <- all[[name]]
    plots[[name]] <- ggplot(D2, aes(x=fADO, y=fEAL, colour = FDR2)) + geom_point(data=D2[D2$FP2 > 0,],aes(size=FP2)) + geom_point(data=D2[D2$FP2 > 0,],aes(size=FP2),shape=1,colour="black") + ggtitle(sprintf("%s - FDR vs FP",name)) + scale + nogrid + bx + by +  scale_size(limits = c(0,306),range=c(0,15))
}

grid.arrange(grobs = plots, ncol=2)



values = seq(0,0.9,by=0.1)
cols = rainbow(10)



###############################################
# plot Spec vs FP
###############################################

plots <- list()
for (name in names){
    D2 <- all[[name]]
    plots[[name]] <- ggplot(D2, aes(x=fADO, y=fEAL,colour = Spec2)) + geom_point(data=D2[D2$FP2 > 0,],aes(size=FP2)) + geom_point(data=D2[D2$FP2 > 0,],aes(size=FP2),shape=1,colour="black") + ggtitle(sprintf("%s - Specificity vs FP",name)) + scale + nogrid + bx + by +  scale_size(limits = c(0,306),range=c(0,15))

}

grid.arrange(grobs = plots, ncol=2)


plots <- list()
for (name in names){
    D2 <- all[[name]]
    plots[[name]] <- ggplot(D2, aes(x=fADO, y=fEAL,colour = Spec2)) + geom_point(data=D2[D2$TN2 > 0,],aes(size=TN2)) + geom_point(data=D2[D2$TN2 > 0,],aes(size=TN2),shape=1,colour="black") + ggtitle(sprintf("%s - Specificity vs TN",name)) + scale + nogrid + bx + by +  scale_size(limits = c(0,306),range=c(0,15))

}

grid.arrange(grobs = plots, ncol=2)





###############################################
# plot sens & spec
###############################################
lw <- 1



par(mfrow=c(2,2), xpd=NA, oma=c(1,0,0,5), cex=0.7)
plot(1, type="n", xlab="fADO", ylab="Sensitivity", xlim=c(0, 0.9), ylim=c(0, 1), main="Sensitivity")

for (name in names){
    D <- t(all[[name]])
    eal0.1 <- which(D["fEAL",]==0.1)
    eal0.5 <- which(D["fEAL",]==0.5)        
    eal0.9 <- which(D["fEAL",]==0.9)    
    lines(as.matrix(D["fADO",eal0]),as.matrix(D["Sens2",eal0.1]),col=coldef2[name],lw=lw)
    lines(as.matrix(D["fADO",eal0.9]),as.matrix(D["Sens2",eal0.9]), lty=3,col=coldef2[name],lw=lw)
    lines(as.matrix(D["fADO",eal0.5]),as.matrix(D["Sens2",eal0.5]), lty=2,col=coldef2[name],lw=lw)    
}
#legend(0.91,1,bty='n', names(coldef2), col=coldef2, lty=1, cex=0.8,lw=1)
#legend(0.91,0.5,bty='n', c("fEAL=0.1","fEAL=0.9"), col="black", lty=c(1,2), cex=0.8,lw=1)


plot(1, type="n", xlab="fADO", ylab="Specificity", xlim=c(0, 0.9), ylim=c(0, 1), main="Specificity")

for (name in names){
    D <- t(all[[name]])
    eal0.1 <- which(D["fEAL",]==0.1)
    eal0.5 <- which(D["fEAL",]==0.5)            
    eal0.9 <- which(D["fEAL",]==0.9)    
    lines(as.matrix(D["fADO",eal0]),as.matrix(D["Spec2",eal0.1]),col=coldef2[name], lw=lw)
    lines(as.matrix(D["fADO",eal0.9]),as.matrix(D["Spec2",eal0.9]), lty=3,col=coldef2[name],lw=lw)
    lines(as.matrix(D["fADO",eal0.5]),as.matrix(D["Spec2",eal0.5]), lty=2,col=coldef2[name],lw=lw)    
}

legend(0.95,1,bty='n', names(coldef2), col=coldef2, lty=1, cex=0.8,lw=lw)
#legend(0.95,0.5,bty='n', c("fEAL=0.1","fEAL=0.9"), col="black", lty=c(1,2), cex=0.8)
legend(0.95,0.5,bty='n', c("fEAL=0.1","fEAL=0.5","fEAL=0.9"), col="black", lty=c(1,2,3), cex=0.8,lw=lw)

# possibly with 3 lines, also eal0.5?


###############################################
# plot FDR both for ADO and EAL
###############################################

plot(1, type="n", xlab="fADO", ylab="FDR", xlim=c(0, 0.9), ylim=c(0, 1), main="FDR vs fADO")

for (name in names){
    D <- t(all[[name]])
    eal0.1 <- which(D["fEAL",]==0.1)
    eal0.9 <- which(D["fEAL",]==0.9)
    eal0.5 <- which(D["fEAL",]==0.5)        
    lines(as.matrix(D["fADO",eal0]),as.matrix(D["FDR2",eal0.1]),col=coldef2[name], lw=lw)
    lines(as.matrix(D["fADO",eal0.9]),as.matrix(D["FDR2",eal0.9]), lty=3,col=coldef2[name],lw=lw)
    lines(as.matrix(D["fADO",eal0.5]),as.matrix(D["FDR2",eal0.5]), lty=2,col=coldef2[name],lw=lw)
}

legend(0.95,1,bty='n', names(coldef2), col=coldef2, lty=1, cex=0.8,lw=lw)
legend(0.95,0.5,bty='n', c("fEAL=0.1","fEAL=0.5","fEAL=0.9"), col="black", lty=c(1,2,3), cex=0.8,lw=lw)


plot(1, type="n", xlab="fEAL", ylab="FDR", xlim=c(0, 0.9), ylim=c(0, 1), main="FDR vs fEAL")

for (name in names){
    D <- t(all[[name]])
    ado0.1 <- which(D["fADO",]==0.1)
    ado0.9 <- which(D["fADO",]==0.9)
    ado0.5 <- which(D["fADO",]==0.5)        
    lines(as.matrix(D["fEAL",ado0]),as.matrix(D["FDR2",ado0.1]),col=coldef2[name], lw=lw)
    lines(as.matrix(D["fEAL",ado0.9]),as.matrix(D["FDR2",ado0.9]), lty=3,col=coldef2[name],lw=lw)
    lines(as.matrix(D["fEAL",ado0.5]),as.matrix(D["FDR2",ado0.5]), lty=2,col=coldef2[name],lw=lw)
}

legend(0.95,1,bty='n', names(coldef2), col=coldef2, lty=1, cex=0.8, lw=lw)
legend(0.95,0.5,bty='n', c("fADO=0.1","fADO=0.5","fADO=0.9"), col="black", lty=c(1,2,3), cex=0.8,lw=lw)


dev.off()






#ggplot(D2, aes(x=fADO, y=fEAL, colour = FDR2)) + geom_point(data=D2[D2$Det2 > 0,],aes(size=Det2)) + geom_point(data=D2[D2$Det2 > 0,],aes(size=Det2),shape=1,colour="black")



#ggsave(plot=p,height=6,width=6,dpi=200, filename="~/example.pdf", useDingbats=FALSE)
#pdf("example.pdf", useDingbats=FALSE)
