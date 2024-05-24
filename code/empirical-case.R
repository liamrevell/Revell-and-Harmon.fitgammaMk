## load packages
library(phytools)
library(geiger)
library(foreach)
library(parallel)
library(doParallel)

## read data from file
squamate.data<-read.csv(file="brandley_table.csv",row.names=1)
rownames(squamate.data)<-gsub(" ","_",rownames(squamate.data))
squamate.tree<-read.nexus(file="squamate.tre")
name.check(squamate.tree,squamate.data)

## fore- and hind-digits
fdigits<-setNames(as.factor(round(squamate.data$Fingers)),
  rownames(squamate.data))
fdigits<-fdigits[squamate.tree$tip.label]
hdigits<-setNames(as.factor(round(squamate.data$Toes)),
  rownames(squamate.data))
hdigits<-hdigits[squamate.tree$tip.label]

## set number of cores to use
niter<-10
ncores<-min(c(parallel::detectCores()-2,
  niter))
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

## fit models (fore-digits)

## ordered two-rate model
ordered1<-matrix(c(
  0,1,0,0,0,0,
  2,0,1,0,0,0,
  0,2,0,1,0,0,
  0,0,2,0,1,0,
  0,0,0,2,0,1,
  0,0,0,0,2,0),6,6,byrow=TRUE)
## standard Mk
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate.tree,fdigits,
    model=ordered1,lik.func="pruning",pi="fitzjohn",
    logscale=TRUE,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
mk_ordered1<-fits[[which(logL==max(logL)[1])]]
plot(mk_ordered1,width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## gamma model
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitgammaMk(squamate.tree,fdigits,
    model=ordered1,pi="fitzjohn",min.alpha=0.001,
    logscale=TRUE,nrates=8,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
gamma_ordered1<-fits[[which(logL==max(logL))[1]]]
plot(as.Qmatrix(gamma_ordered1),width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)

## ordered ten-rate model
ordered2<-matrix(c(
  0,1,0,0,0,0,
  2,0,3,0,0,0,
  0,4,0,5,0,0,
  0,0,6,0,7,0,
  0,0,0,8,0,9,
  0,0,0,0,10,0),6,6,byrow=TRUE)
## standard Mk
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate.tree,fdigits,
    model=ordered2,lik.func="pruning",pi="fitzjohn",
    logscale=TRUE,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
mk_ordered2<-fits[[which(logL==max(logL))[1]]]
plot(mk_ordered2,width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## gamma model
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitgammaMk(squamate.tree,fdigits,
    model=ordered2,pi="fitzjohn",min.alpha=0.001,
    logscale=TRUE,nrates=8,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
gamma_ordered2<-fits[[which(logL==max(logL))[1]]]
plot(as.Qmatrix(gamma_ordered2),width=TRUE,color=TRUE)

## directional one-rate model
directional1<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,1,0,0,0,0,
  0,0,1,0,0,0,
  0,0,0,1,0,0,
  0,0,0,0,1,0),6,6,byrow=TRUE)
## standard Mk
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate.tree,fdigits,
    model=directional1,lik.func="pruning",pi="fitzjohn",
    logscale=TRUE,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
mk_directional1<-fits[[which(logL==max(logL)[1])]]
plot(mk_directional1,width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## gamma model
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitgammaMk(squamate.tree,fdigits,
    model=directional1,pi="fitzjohn",min.alpha=0.001,
    logscale=TRUE,nrates=8,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
gamma_directional1<-fits[[which(logL==max(logL))[1]]]
plot(as.Qmatrix(gamma_directional1),width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)

## direction five-rate model
directional2<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,2,0,0,0,0,
  0,0,3,0,0,0,
  0,0,0,4,0,0,
  0,0,0,0,5,0),6,6,byrow=TRUE)
## standard Mk
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate.tree,fdigits,
    model=directional2,lik.func="pruning",pi="fitzjohn",
    logscale=TRUE,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
mk_directional2<-fits[[which(logL==max(logL)[1])]]
plot(mk_directional2,width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## gamma model
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitgammaMk(squamate.tree,fdigits,
    model=directional2,pi="fitzjohn",min.alpha=0.001,
    logscale=TRUE,nrates=8,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
gamma_directional2<-fits[[which(logL==max(logL))[1]]]
plot(as.Qmatrix(gamma_directional2),width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)

stopCluster(cl=mc)

marginal.gamma_ordered2<-fitgammaMk(squamate.tree,
  fdigits,fixedQ=as.Qmatrix(gamma_ordered2),
  alpha.init=gamma_ordered2$alpha,
  pi="fitzjohn",marginal=TRUE)


save(mk_directional1,gamma_directional1,
  mk_ordered1,gamma_ordered1,
  mk_directional2,gamma_directional2,
  mk_ordered2,gamma_ordered2,
  marginal.gamma_ordered2,
  file="squamate-fdigit-models.rda")

anova(mk_directional1,gamma_directional1,
  mk_ordered1,gamma_ordered1,
  mk_directional2,gamma_directional2,
  mk_ordered2,gamma_ordered2)

ace.gamma_ordered2<-ancr(gamma_ordered2)

cols<-setNames(
  viridisLite::viridis(n=ncol(ace.gamma_ordered2$ace)),
  colnames(ace.gamma_ordered2$ace))
plot(marginal.gamma_ordered2,
  colors=c("black","navy","blue","lightblue","ivory2"),
  ftype="i",offset=0.5)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
for(i in 1:Ntip(squamate.tree)){
  polygon(x=pp$xx[i]+c(0,3,3,0),
    y=pp$yy[i]+c(-0.5,-0.5,0.5,0.5),
    col=cols[as.numeric(fdigits[i])],border=FALSE)
}
pie_cex<-apply(ace.gamma_ordered2$ace,1,
  function(x) if(any(x>0.95)) 0.3 else 0.5)
par(fg="transparent")
nodelabels(pie=ace.gamma_ordered2$ace,piecol=cols,
  cex=pie_cex)
par(fg="black")
legend(x=mean(par()$usr[1:2]),y=5,
  colnames(ace.gamma_ordered2$ace),
  pch=16,pt.cex=2,col=cols,horiz=TRUE,
  bty="n",xjust=0.5)


## set number of cores to use
niter<-10
ncores<-min(c(parallel::detectCores()-2,
  niter))
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

## fit models (hind-digits)
## ordered two-rate model
ordered1<-matrix(c(
  0,1,0,0,0,0,
  2,0,1,0,0,0,
  0,2,0,1,0,0,
  0,0,2,0,1,0,
  0,0,0,2,0,1,
  0,0,0,0,2,0),6,6,byrow=TRUE)
## standard Mk
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate.tree,hdigits,
    model=ordered1,lik.func="pruning",pi="fitzjohn",
    logscale=TRUE,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
mk_ordered1.h<-fits[[which(logL==max(logL))[1]]]
plot(mk_ordered1.h,width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## gamma model
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitgammaMk(squamate.tree,hdigits,
    model=ordered1,pi="fitzjohn",min.alpha=0.001,
    logscale=TRUE,nrates=8,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
gamma_ordered1.h<-fits[[which(logL==max(logL))[1]]]
plot(as.Qmatrix(gamma_ordered1.h),width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## ordered ten-rate model
ordered2<-matrix(c(
  0,1,0,0,0,0,
  2,0,3,0,0,0,
  0,4,0,5,0,0,
  0,0,6,0,7,0,
  0,0,0,8,0,9,
  0,0,0,0,10,0),6,6,byrow=TRUE)
## standard Mk
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate.tree,hdigits,
    model=ordered2,lik.func="pruning",pi="fitzjohn",
    logscale=TRUE,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
mk_ordered2.h<-fits[[which(logL==max(logL))[1]]]
plot(mk_ordered2.h,width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## gamma model
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitgammaMk(squamate.tree,hdigits,
    model=ordered2,pi="fitzjohn",min.alpha=0.001,
    logscale=TRUE,nrates=8,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
gamma_ordered2.h<-fits[[which(logL==max(logL))[1]]]
plot(as.Qmatrix(gamma_ordered2.h),width=TRUE,color=TRUE)


## directional one-rate model
directional1<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,1,0,0,0,0,
  0,0,1,0,0,0,
  0,0,0,1,0,0,
  0,0,0,0,1,0),6,6,byrow=TRUE)
## standard Mk
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate.tree,hdigits,
    model=directional1,lik.func="pruning",pi="fitzjohn",
    logscale=TRUE,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
mk_directional1.h<-fits[[which(logL==max(logL))[1]]]
plot(mk_directional1.h,width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## gamma model
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitgammaMk(squamate.tree,hdigits,
    model=directional1,pi="fitzjohn",min.alpha=0.001,
    logscale=TRUE,nrates=8,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
gamma_directional1.h<-fits[[which(logL==max(logL))[1]]]
plot(as.Qmatrix(gamma_directional1.h),width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## direction five-rate model
directional2<-matrix(c(
  0,0,0,0,0,0,
  1,0,0,0,0,0,
  0,2,0,0,0,0,
  0,0,3,0,0,0,
  0,0,0,4,0,0,
  0,0,0,0,5,0),6,6,byrow=TRUE)
## standard Mk
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitMk(squamate.tree,hdigits,
    model=directional2,lik.func="pruning",pi="fitzjohn",
    logscale=TRUE,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
mk_directional2.h<-fits[[which(logL==max(logL)[1])]]
plot(mk_directional2.h,width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
## gamma model
fits<-foreach(i=1:niter)%dopar%{
  phytools::fitgammaMk(squamate.tree,hdigits,
    model=directional2,pi="fitzjohn",min.alpha=0.001,
    logscale=TRUE,nrates=8,rand_start=TRUE)
}
logL<-sapply(fits,logLik)
logL
gamma_directional2.h<-fits[[which(logL==max(logL))[1]]]
plot(as.Qmatrix(gamma_directional2.h),width=TRUE,color=TRUE,
  mar=rep(0.1,4),offset=0.025)
stopCluster(cl=mc)

marginal.gamma_ordered2.h<-fitgammaMk(squamate.tree,
  fdigits,fixedQ=as.Qmatrix(gamma_ordered2.h),
  alpha.init=gamma_ordered2.h$alpha,
  pi="fitzjohn",marginal=TRUE)

save(mk_directional1.h,gamma_directional1.h,
  mk_ordered1.h,gamma_ordered1.h,
  mk_directional2.h,gamma_directional2.h,
  mk_ordered2.h,gamma_ordered2.h,
  marginal.gamma_ordered2.h,
  file="squamate-hdigit-models.rda")

aov_table.h<-anova(mk_directional1.h,gamma_directional1.h,
  mk_ordered1.h,gamma_ordered1.h,
  mk_directional2.h,gamma_directional2.h,
  mk_ordered2.h,gamma_ordered2.h)



ace.gamma_ordered2.h<-ancr(gamma_ordered2.h)
cols<-setNames(
  viridisLite::viridis(n=ncol(ace.gamma_ordered2.h$ace)),
  colnames(ace.gamma_ordered2.h$ace))
plot(marginal.gamma_ordered2.h,
  colors=c("black","navy","blue","lightblue","ivory2"),
  ftype="i",offset=0.4,lwd=1,fsize=0.2)
pp<-get("last_plot.phylo",envir=.PlotPhyloEnv)
for(i in 1:Ntip(squamate.tree)){
  polygon(x=pp$xx[i]+c(0,3,3,0),
    y=pp$yy[i]+c(-0.5,-0.5,0.5,0.5),
    col=cols[as.numeric(hdigits[i])],border=FALSE)
}
pie_cex<-apply(ace.gamma_ordered2.h$ace,1,
  function(x) if(any(x>0.95)) 0.3 else 0.5)
par(fg="transparent")
nodelabels(pie=ace.gamma_ordered2.h$ace,piecol=cols,
  cex=pie_cex)
par(fg="black")
legend(x=mean(par()$usr[1:2]),y=5,
  colnames(ace.gamma_ordered2.h$ace),
  pch=16,pt.cex=2,col=cols,horiz=TRUE,
  bty="n",xjust=0.5)
