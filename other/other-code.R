library(phytools)
tree<-pbtree(n=400)
x<-fastBM(tree,sig2=0.5,mu=40)
y<-as.factor(round(x))
Y<-to.matrix(y,seq(min(round(x)),max(round(x))))
MODEL<-matrix(0,ncol(Y),ncol(Y),dimnames=list(colnames(Y),colnames(Y)))
for(i in 1:ncol(MODEL)){
  if(i<nrow(MODEL)) MODEL[i,i+1]<-1
  if(i>1) MODEL[i,i-1]<-2
}
MODEL

fitMk(tree,Y,model=MODEL,pi="fitzjohn")

a<-setNames(c(0,0,0),LETTERS[1:3])
b<-setNames(c(0,0,1),LETTERS[1:3])
c<-setNames(c(0,1,0),LETTERS[1:3])
d<-setNames(c(1,0,0),LETTERS[1:3])
e<-setNames(c(0,1,1),LETTERS[1:3])
f<-setNames(c(1,0,1),LETTERS[1:3])
g<-setNames(c(1,1,0),LETTERS[1:3])
h<-setNames(c(1,1,1),LETTERS[1:3])

exp(logLik(fitMk(tree,to.matrix(a,0:1),fixedQ=Q)))+
  exp(logLik(fitMk(tree,b,fixedQ=Q)))+
  exp(logLik(fitMk(tree,c,fixedQ=Q)))+
  exp(logLik(fitMk(tree,d,fixedQ=Q)))+
  exp(logLik(fitMk(tree,e,fixedQ=Q)))+
  exp(logLik(fitMk(tree,f,fixedQ=Q)))+
  exp(logLik(fitMk(tree,g,fixedQ=Q)))+
  exp(logLik(fitMk(tree,to.matrix(h,0:1),fixedQ=Q)))
  

mat<-matrix(c(1:10,rep(11,5)),3,5,byrow=TRUE)
set.seed(17)
smap<-make.simmap(tree,x,Q=Q,nsim=10,message=FALSE)
layout(mat,heights=c(0.45,0.45,0.1))
cols<-setNames(c("white","black"),0:1)
plot(smap,cols,ftype="i",lwd=5,outline=TRUE)
par(mar=rep(0,4))
plot(NA,xlim=c(0,1),ylim=c(0,1),bty="n",axes=FALSE,
  xlab="",ylab="")
legend("center",legend=0:1,pch=22,pt.bg=c("white","black"),
  horiz=TRUE,bty="n",cex=1.5,pt.cex=2)
