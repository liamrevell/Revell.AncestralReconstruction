library(phytools)
library(doParallel)
library(foreach)

set.seed(99)
ntaxa<-1001
nsim<-100
trees<-pbtree(n=ntaxa,scale=1,nsim=nsim)

alpha<-1
nrates<-10
r<-qgamma(seq(1/(2*nrates),1,by=1/nrates),alpha,alpha)
r<-r/mean(r)

sim_trees<-trees
for(i in 1:nsim){
  sim_trees[[i]]$edge.length<-sim_trees[[i]]$edge.length*
    sample(r,nrow(sim_trees[[i]]$edge),replace=TRUE)
}

Q<-model<-matrix(c(0,1,0,1,0,1,0,1,0),3,3,
  dimnames=list(0:2,0:2))
#Q<-model<-matrix(c(0,1,1,0),2,2,
#  dimnames=list(0:1,0:1))
diag(Q)<--rowSums(Q)

y<-lapply(sim_trees,sim.Mk,Q=Q,internal=TRUE)
y.obs<-lapply(y,function(x,ntaxa) x[1:ntaxa],ntaxa=ntaxa)

ncores<-min(nsim,detectCores()-1)
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

mk_fits<-foreach(i=1:nsim)%dopar%{
  phytools::fitMk(tree=trees[[i]],x=y.obs[[i]],
    model=model,pi="fitzjohn",lik.func="pruning",
    rand_start=TRUE)
}

# mk_fits<-mapply(fitMk,tree=trees,x=y.obs,
#  MoreArgs=list(model=model,pi="fitzjohn",
#    lik.func="pruning"),SIMPLIFY=FALSE)

gamma_fits<-foreach(i=1:nsim)%dopar%{
  phytools::fitgammaMk(tree=trees[[i]],x=y.obs[[i]],
    model=model,pi="fitzjohn",rand_start=TRUE,nrates=10)
}

# gamma_fits<-mapply(fitgammaMk,tree=trees,x=y.obs,
#  MoreArgs=list(model=model,pi="fitzjohn"),
#  SIMPLIFY=FALSE)

mapply(anova,mk_fits,gamma_fits,SIMPLIFY=FALSE)

exp(mean(log(sapply(gamma_fits,function(x) x$alpha))))

Alpha<-sapply(gamma_fits,function(x) x$alpha)

hist(log(Alpha),breaks=40)
abline(v=log(0.5),lwd=3)

a<-lapply(y,function(x,ntaxa) x[-c(1:ntaxa)],ntaxa=ntaxa)

anc_y<-lapply(mk_fits,ancr)
anc_y<-foreach(i=1:nsim)%dopar%{ phytools::ancr(mk_fits[[i]])}
asr_y<-data.frame(
  bin=seq(0.01,0.99,by=0.02),
  count=rep(0,50),
  actual=rep(0,50),
  prop=rep(0,50))
for(i in 1:nsim){
  for(j in 1:nrow(asr_y)){
    ii<-intersect(which(anc_y[[i]]$ace[,2]>=(asr_y[j,1]-0.01)),
      which(anc_y[[i]]$ace[,2]<=(asr_y[j,1]+0.01)))
    asr_x[j,"count"]<-asr_y[j,"count"]+length(ii)
    asr_x[j,"actual"]<-asr_y[j,"actual"]+sum(x.obs[[i]][ntaxa+ii]==1)
  }
}
asr_x[,"prop"]<-asr_x[,"actual"]/asr_x[,"count"]
asr_x[is.nan(asr_x[,"prop"]),"prop"]<-0



stopCluster(mc) ## stop cluster
