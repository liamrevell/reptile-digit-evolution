library(phytools)
library(geiger)
squamate.data<-read.csv(file="brandley_table.csv",row.names=1)
rownames(squamate.data)<-gsub(" ","_",rownames(squamate.data))
squamate.tree<-read.nexus(file="squamate.tre")
name.check(squamate.tree,squamate.data)

ff<-floor(squamate.data$Toes)
cc<-ceiling(squamate.data$Toes)

foo<-function(x,y){
    obj<-intersect(x,y)
    if(length(obj)==0) obj<-paste(x,y,sep="+")
    obj
}

toes<-mapply(foo,x=ff,y=cc)

toes<-setNames(as.factor(ff),rownames(squamate.data))

toes_mk<-fitHRM(squamate.tree,toes,umbral=TRUE,ncat=1,niter=20,
  parallel=TRUE,ncores=10,pi="fitzjohn",ordered=TRUE,order=levels(toes),
  lik.func="pruning",logscale=TRUE)
plot(as.Qmatrix(toes_mk),show.zeros=FALSE)

toes_hrm<-fitHRM(squamate.tree,toes,umbral=TRUE,ncat=c(2,2,2,2,2,2),niter=10,
  parallel=TRUE,ncores=10,pi="fitzjohn",ordered=TRUE,order=levels(toes),
  lik.func="pruning",logscale=TRUE)
sapply(toes_hrm$all.fits,logLik)4

getwd()
setwd("../courses/Biol634-spring2023/project/")
list.files()
setwd("doi_10.5061_dryad.4d5h8g12__v1/")
list.files()

squamate.tree<-read.tree(file="Species_level_supertree.tre")
squamate.data<-read.csv(file="Morphometrics.csv",row.names=1)
chk<-name.check(squamate.tree,squamate.data)

squamate.tree<-drop.tip(squamate.tree,chk$tree_not_data)
squamate.tree

name.check(squamate.tree,squamate.data)

hind_digits<-setNames(as.factor(squamate.data$HDp1-1),rownames(squamate.data))

toes_mk<-fitHRM(squamate.tree,hind_digits,umbral=TRUE,ncat=1,niter=10,
  parallel=TRUE,ncores=10,pi="fitzjohn",ordered=TRUE,order=levels(hind_digits),
  lik.func="pruning",logscale=TRUE)

options(scipen=6)
plot(toes_mk,signif=4)

plot(as.Qmatrix(toes_mk),show.zeros=FALSE,signif=4,
  width=TRUE,color=TRUE)

toes_hrm<-fitHRM(squamate.tree,hind_digits,umbral=TRUE,ncat=c(2,1,1,1,1,2),niter=10,
  parallel=TRUE,ncores=10,pi="fitzjohn",ordered=TRUE,order=levels(hind_digits),
  lik.func="pruning",logscale=TRUE)

plot(toes_hrm)

anova(toes_mk,toes_hrm)


plotTree(squamate.tree)

test<-read.tree(file="../../../../manuscripts/Revell-and-Mahler.reptile-viviparity/ele12168-sup-0006-data_file_1.txt")
