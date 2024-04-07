## load packages
library(doParallel)
library(foreach)

set.seed(99)
ntaxa<-501
nsim<-100

## make cluster
ncores<-min(nsim,detectCores()-1)
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

y<-lapply(trees,fastBM,bounds=c(-4,4),internal=TRUE)
y.obs<-lapply(x,function(x,ntaxa) x[1:ntaxa],ntaxa=ntaxa)

anc_bm<-foreach(i=1:nsim)%dopar%{
  phytools::fastAnc(tree=trees[[i]],x=y.obs[[i]],CI=TRUE)
}

bounded_fits<-mapply(bounded_bm,tree=trees,x=y.obs,
  lims=lapply(y.obs,range),
  MoreArgs=list(lik.func="eigen",parallel=TRUE),
  SIMPLIFY=FALSE)
anc_bounded<-foreach(i=1:nsim)%dopar%{
  phytools::ancr(bounded_fits[[i]])
}

a<-lapply(y,function(x,ntaxa) x[-c(1:ntaxa)],ntaxa=ntaxa)

foo<-function(obj,a){
  sum(mapply('&&',a>=obj$CI95[,1],a<=obj$CI95[,2]))/length(a)
}
onCI.bm<-mapply(foo,obj=anc_bm,a=a,SIMPLIFY=TRUE)
onCI.bm
onCI.bounded<-mapply(foo,obj=anc_bounded,a=a,SIMPLIFY=TRUE)
onCI.bounded

stopCluster(mc)

save(trees,y,y.obs,a,anc_bm,anc_bounded,
  bounded_fits,nsim,ntaxa,onCI.bm,onCI.bounded,
  file="asr-bounded.rda")

