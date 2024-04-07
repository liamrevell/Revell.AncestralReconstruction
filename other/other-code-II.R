```{r, echo=FALSE, fig.width=7, fig.height=5, dpi=300, fig.cap="See main text for more details.", out.width = "100%"}
foo<-function(n) gray.colors(n,end=0,start=1)
par(cex.axis=0.8)
filled.contour(x=x0,y=x1,z=L,nlevels=50,las=1,bty="n",
  cex.axis=0.8,xlab="root node",ylab="internal node",
  color.palette=foo)
abline(v=fit$par[1])
arrows(x0=fit$par[1],y0=fit$par[2],x1=fit$par[1],y1=0,lwd=1,col=grey(0.5),lend=3,length=0.1)
arrows(x0=fit$par[1],y0=fit$par[2],x1=par()$usr[1],y1=fit$par[2],lwd=1,col=grey(0.5),lend=3,length=0.1)
```
```{r, echo=FALSE, fig.width=7, fig.height=5, dpi=300, fig.cap="See main text for more details.", out.width = "100%"}
image(x0,x1,L,col=grey.colors(n=1000,1,0))
```