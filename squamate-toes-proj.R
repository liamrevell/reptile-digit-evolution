library(phytools)
library(geiger)

squamate.tree<-read.tree(file="Species_level_supertree.tre")
squamate.data<-read.csv(file="Morphometrics.csv",row.names=1)
chk<-name.check(squamate.tree,squamate.data)

squamate.tree<-drop.tip(squamate.tree,chk$tree_not_data)
squamate.tree

name.check(squamate.tree,squamate.data)

## fitting an ordered and hidden-rates model to hind-digit number

hind_digits<-setNames(as.factor(squamate.data$HDp1-1),rownames(squamate.data))

toes_mk<-fitHRM(squamate.tree,hind_digits,umbral=TRUE,ncat=1,niter=10,
  parallel=TRUE,ncores=10,pi="fitzjohn",ordered=TRUE,order=levels(hind_digits),
  lik.func="pruning",logscale=TRUE)

png(file="mk_rear_toes.png",width=10,height=6,units="in",res=300)
options(scipen=6)
plot(toes_mk,signif=5,ylim=c(-1,1),asp=0.4)
dev.off()

png(file="mk_rear_toes-2.png",width=12,height=8,units="in",res=300)
plot(as.Qmatrix(toes_mk),show.zeros=FALSE,signif=5,
  width=TRUE,color=TRUE,xlim=c(-1.8,1))
dev.off()

toes_hrm<-fitHRM(squamate.tree,hind_digits,umbral=TRUE,ncat=c(2,1,1,1,1,2),niter=10,
  parallel=TRUE,ncores=10,pi="fitzjohn",ordered=TRUE,order=levels(hind_digits),
  lik.func="pruning",logscale=TRUE)

png(file="hrm_rear_toes.png",width=10,height=6,units="in",res=300)
plot(toes_hrm,asp=0.5,signif=5)
dev.off()

anova(toes_mk,toes_hrm)

toes_mk_asr<-ancr(toes_mk)

png(file="mk_toes_asr.png",width=12,height=8,units="in",res=300)
h<-max(nodeHeights(squamate.tree))
plotTree(squamate.tree,direction="upwards",ftype="off",lwd=1,
  color="black",ylim=c(0,1.05*h))
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cols<-setNames(viridisLite::viridis(n=6,direction=-1),
  levels(hind_digits))
for(i in 1:Ntip(squamate.tree)){
  COL<-cols[hind_digits[squamate.tree$tip.label[i]]]
  polygon(pp$xx[i]+c(-0.5,-0.5,0.5,0.5),
    pp$yy[i]+c(0,0.05*h,0.05*h,0),
    col=COL,border=FALSE)
}
legend("bottomleft",levels(hind_digits),pch=22,pt.bg=cols)
par(fg="transparent")
cex_pies<-apply(toes_mk_asr$ace,1,function(x) if(max(x>0.9)) 0.2 else 0.4)
nodelabels(pie=toes_mk_asr$ace,piecol=cols,cex=cex_pies)
par(fg="black")
dev.off()

toes_hrm_asr<-ancr(toes_hrm)

png(file="hrm_toes_asr.png",width=12,height=8,units="in",res=300)
h<-max(nodeHeights(squamate.tree))
plotTree(squamate.tree,direction="upwards",ftype="off",lwd=1,
         color="black",ylim=c(0,1.05*h))
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cols<-setNames(viridisLite::viridis(n=6,direction=-1),
               levels(hind_digits))
for(i in 1:Ntip(squamate.tree)){
  COL<-cols[hind_digits[squamate.tree$tip.label[i]]]
  polygon(pp$xx[i]+c(-0.5,-0.5,0.5,0.5),
          pp$yy[i]+c(0,0.05*h,0.05*h,0),
          col=COL,border=FALSE)
}
legend("bottomleft",levels(hind_digits),pch=22,pt.bg=cols)
par(fg="transparent")
cex_pies<-apply(hide.hidden(toes_hrm_asr),1,function(x) if(max(x>0.9)) 0.2 else 0.4)
nodelabels(pie=hide.hidden(toes_hrm_asr),piecol=cols,cex=cex_pies)
par(fg="black")
dev.off()

## fitting an ordered and hidden-rates model to hind-digit number

fore_digits<-setNames(as.factor(squamate.data$FDp1-1),rownames(squamate.data))

fingers_mk<-fitHRM(squamate.tree,fore_digits,umbral=TRUE,ncat=1,niter=10,
  parallel=TRUE,ncores=10,pi="fitzjohn",ordered=TRUE,order=levels(fore_digits),
  lik.func="pruning",logscale=TRUE)

options(scipen=6)
plot(fingers_mk,signif=4)

plot(as.Qmatrix(fingers_mk),show.zeros=FALSE,signif=4,
  width=TRUE,color=TRUE)

fingers_hrm<-fitHRM(squamate.tree,fore_digits,umbral=TRUE,ncat=c(2,1,1,1,1,2),niter=10,
  parallel=TRUE,ncores=10,pi="fitzjohn",ordered=TRUE,order=levels(fore_digits),
  lik.func="pruning",logscale=TRUE)

plot(fingers_hrm,signif=6)

anova(fingers_mk,fingers_hrm)

fingers_mk_asr<-ancr(fingers_mk)

png(file="mk_fingers_asr.png",width=12,height=8,units="in",res=600)
h<-max(nodeHeights(squamate.tree))
plotTree(squamate.tree,direction="upwards",ftype="off",lwd=1,
         color="black",ylim=c(0,1.05*h))
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cols<-setNames(viridisLite::viridis(n=6,direction=-1),
               levels(hind_digits))
for(i in 1:Ntip(squamate.tree)){
  COL<-cols[hind_digits[squamate.tree$tip.label[i]]]
  polygon(pp$xx[i]+c(-0.5,-0.5,0.5,0.5),
          pp$yy[i]+c(0,0.05*h,0.05*h,0),
          col=COL,border=FALSE)
}
legend("bottomleft",levels(hind_digits),pch=22,pt.bg=cols)
par(fg="transparent")
cex_pies<-apply(fingers_mk_asr$ace,1,function(x) if(max(x>0.9)) 0.2 else 0.4)
nodelabels(pie=fingers_mk_asr$ace,piecol=cols,cex=cex_pies)
par(fg="black")
dev.off()

fingers_hrm_asr<-ancr(fingers_hrm)

png(file="hrm_fingers_asr.png",width=12,height=8,units="in",res=300)
h<-max(nodeHeights(squamate.tree))
plotTree(squamate.tree,direction="upwards",ftype="off",lwd=1,
         color="black",ylim=c(0,1.05*h))
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
cols<-setNames(viridisLite::viridis(n=6,direction=-1),
               levels(hind_digits))
for(i in 1:Ntip(squamate.tree)){
  COL<-cols[hind_digits[squamate.tree$tip.label[i]]]
  polygon(pp$xx[i]+c(-0.5,-0.5,0.5,0.5),
          pp$yy[i]+c(0,0.05*h,0.05*h,0),
          col=COL,border=FALSE)
}
legend("bottomleft",levels(hind_digits),pch=22,pt.bg=cols)
par(fg="transparent")
cex_pies<-apply(hide.hidden(fingers_hrm_asr),1,function(x) if(max(x>0.9)) 0.2 else 0.4)
nodelabels(pie=hide.hidden(fingers_hrm_asr),piecol=cols,cex=cex_pies)
par(fg="black")
dev.off()
