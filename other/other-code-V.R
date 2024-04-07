## load packages
library(doParallel)
library(foreach)

## part 1

set.seed(99)
ntaxa<-501
nsim<-100

## make cluster
ncores<-min(nsim,detectCores()-1)
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

Q<-matrix(c(
  -0.2, 0.2, 0.0,
  0.2,-0.3, 0.1,
  0.0, 0.0, 0.0),3,3,byrow=TRUE,
  dimnames=list(c("0","1","1*"),c("0","1","1*")))
foo<-function(tree,Q){
  x<-sim.Mk(tree,Q,internal=TRUE,anc="0")
  while(("0"%in%as.character(x[tree$tip.label]))==FALSE)
    x<-sim.Mk(tree,Q,internal=TRUE,anc="0")
  x
}
x<-lapply(trees,foo,Q=Q)
foo<-function(x){
  y<-as.character(x)
  y[y=="1*"]<-"1"
  setNames(as.factor(y),names(x))
}
x.obs<-lapply(x,foo)
fits_x<-foreach(i=1:nsim)%dopar%{
  phytools::fitMk(trees[[i]],x[[i]][1:ntaxa],model="ARD")
}
# fits_x<-mapply(function(t,x) fitMk(t,x[t$tip.label],
#   model="ARD"),t=trees,x=x.obs,SIMPLIFY=FALSE)
anc_x<-lapply(fits_x,ancr)
asr_x<-data.frame(
  bin=seq(0.01,0.99,by=0.02),
  count=rep(0,50),
  actual=rep(0,50),
  prop=rep(0,50))
for(i in 1:nsim){
  for(j in 1:nrow(asr_x)){
    ii<-intersect(which(anc_x[[i]]$ace[,2]>=(asr_x[j,1]-0.01)),
      which(anc_x[[i]]$ace[,2]<=(asr_x[j,1]+0.01)))
    asr_x[j,"count"]<-asr_x[j,"count"]+length(ii)
    asr_x[j,"actual"]<-asr_x[j,"actual"]+sum(x.obs[[i]][ntaxa+ii]==1)
  }
}
asr_x[,"prop"]<-asr_x[,"actual"]/asr_x[,"count"]
asr_x[is.nan(asr_x[,"prop"]),"prop"]<-0

stopCluster(mc)

## part 2


fits_x.hrm<-mapply(function(t,x) fitHRM(t,x[t$tip.label],
  ncat=c(1,2),parallel=TRUE),t=trees,x=x.obs,SIMPLIFY=FALSE)
anc_x.hrm<-lapply(fits_x.hrm,ancr)
anc_x.hidden<-lapply(anc_x.hrm,hide.hidden)
asr_x.hrm<-data.frame(
  bin=seq(0.01,0.99,by=0.02),
  count=rep(0,50),
  actual=rep(0,50),
  prop=rep(0,50))
for(i in 1:nsim){
  for(j in 1:nrow(asr_x.hrm)){
    ii<-intersect(which(anc_x.hidden[[i]][,2]>=(asr_x.hrm[j,1]-0.01)),
      which(anc_x.hidden[[i]][,2]<=(asr_x.hrm[j,1]+0.01)))
    asr_x.hrm[j,"count"]<-asr_x.hrm[j,"count"]+length(ii)
    asr_x.hrm[j,"actual"]<-asr_x.hrm[j,"actual"]+sum(x.obs[[i]][ntaxa+ii]==1)
  }
}
asr_x.hrm[,"prop"]<-asr_x.hrm[,"actual"]/asr_x.hrm[,"count"]
asr_x.hrm[is.nan(asr_x.hrm[,"prop"]),"prop"]<-0