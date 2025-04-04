par(mfrow=c(1,2))

## set seed
set.seed(99)
## load packages
library(foreach)
library(parallel)
library(doParallel)
## LR test function
lr.test<-function(obj1,obj2){
  l1<-logLik(obj1)
  l2<-logLik(obj2)
  LR<-2*(l2-l1)
  as.numeric(
    pchisq(LR,df=attr(l2,"df")-attr(l1,"df"),
      lower.tail=FALSE))
}
## simulation conditions
niter<-3
nrep<-100 # 100
alpha<-c(0.25,0.5,1,2,4,8)
ntaxa<-2000 # 2000
## set number of cores to use
ncores<-min(c(parallel::detectCores()-2,
  niter))
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)
## results matrix
alphaRESULTS<-matrix(NA,nrep*length(alpha),6,
  dimnames=list(1:(nrep*length(alpha)),
    c("alpha","est.alpha","est.q","logL",
      "logL(null)","P")))
## generate random tree(s)
trees<-pbtree(n=ntaxa,scale=10,nsim=nrep)
## generate Q matrix for simulation
Q<-matrix(c(
  -0.5, 0.5, 0.0,
  0.5,-1.0, 0.5,
  0.0, 0.5,-0.5),3,3,
  dimnames=list(letters[1:3],letters[1:3]))
## design matrix of model
MODEL<- matrix(c(
  0,1,0,
  1,0,1,
  0,1,0),3,3,
  dimnames=list(letters[1:3],letters[1:3]))
## lists to store generating conditions/data
rates<-list()
x<-list()
for(i in 1:length(alpha)){
  rates[[i]]<-list()
  x[[i]]<-list()
  for(j in 1:nrep){
    max.alpha<-1000
    cat("\n#################\n\n")
    ## simulate data
    rates[[i]][[j]]<-rgamma(n=nrow(trees[[j]]$edge),alpha[i],alpha[i])
    sim_tree<-trees[[j]]
    sim_tree$edge.length<-trees[[j]]$edge.length*rates[[i]][[j]]
    x[[i]][[j]]<-sim.Mk(sim_tree,Q)
    X<-to.matrix(x[[i]][[j]],letters[1:ncol(Q)])
    plotTree(sim_tree,ftype="off",lwd=1,mar=c(1.1,1.1,4.1,1.1))
    tmp<-bquote(alpha==.(round(alpha[i],3)))
    title(main=tmp)
    fits<-foreach(k=1:niter)%dopar%{
      phytools::fitMk(trees[[j]],X,model=MODEL,
        pi="fitzjohn",rand_start=TRUE)
    }
    logL<-sapply(fits,logLik)
    cat("null fits:\n")
    print(logL)
    fit_null<-fits[[which.max(logL)[1]]]
    fits<-foreach(k=1:niter)%dopar%{
      phytools::fitgammaMk(trees[[j]],X,model=MODEL,min.alpha=0.001,
        pi="fitzjohn",rand_start=TRUE,nrates=8)
    }
    logL<-sapply(fits,logLik)
    cat("gamma fits:\n")
    print(logL)
    best_fit<-fits[[which.max(logL)[1]]]
    while((abs(best_fit$alpha-max.alpha)<1e-08)&(max.alpha<=1e5)){
      max.alpha<-max.alpha*10
      fits<-c(foreach(k=1:niter)%dopar%{
        phytools::fitgammaMk(trees[[j]],X,model=MODEL,min.alpha=0.001,
          pi="fitzjohn",rand_start=TRUE,nrates=8,max.alpha=max.alpha)
      },fits)
      logL<-sapply(fits,logLik)
      cat("gamma fits:\n")
      print(logL)
      best_fit<-fits[[which.max(logL)[1]]]
    }
    counter<-0
    while((logLik(fit_null)>logLik(best_fit))&(counter<3)){
      counter<-counter+1
      fits<-c(foreach(k=1:niter)%dopar%{
        phytools::fitgammaMk(trees[[j]],X,model=MODEL,min.alpha=0.001,
          pi="fitzjohn",rand_start=TRUE,nrates=8,max.alpha=max.alpha)
      },fits)
      logL<-sapply(fits,logLik)
      cat("gamma fits:\n")
      print(logL)
      best_fit<-fits[[which.max(logL)[1]]]
    }
    alphaRESULTS[(i-1)*nrep+j,]<-
      c(alpha[i],best_fit$alpha,best_fit$rates,
        logLik(best_fit),logLik(fit_null),
        lr.test(fit_null,best_fit))
    print(alphaRESULTS[(i-1)*nrep+j,,drop=FALSE])
    cat("\n")
  }
}
stopCluster(cl=mc)

alpha<-c(0.25,0.5,1,2,4,8)
par(bty="n",mfrow=c(1,2),mar=c(5.1,4.1,2.1,1.1))
boxplot(alphaRESULTS[,"est.alpha"]~as.factor(alphaRESULTS[,"alpha"]),log="y",
  bty="n",xlab=expression(paste(alpha," shape parameter")),
  ylab=expression(paste("estimated value of ",alpha)),las=1,cex.axis=0.8,
  xlim=c(-0.5,6.5),ylim=c(0.1,1e6))
abline(h=alpha,lty="dashed",col="blue")
boxplot(alphaRESULTS[,"est.alpha"]~as.factor(alphaRESULTS[,"alpha"]),log="y",
  add=TRUE,cex=1.2,pch=21,bg="grey",las=1,cex.axis=0.8)
for(i in 1:length(alpha)){
  tmp<-bquote(alpha==.(round(alpha[i],3)))
  text(x=-0.7,y=1.2*alpha[i],tmp,pos=4,cex=0.8)
}
legend("topleft",expression(paste("generating value of ",alpha)),
  lwd=1,col="blue",lty="dashed",bty="n",cex=0.8)
mtext("A",line=0.5,adj=-0.1,font=2,cex=1.5)
foo<-function(p) qchisq(p,df=1,lower.tail=FALSE)/2
boxplot(foo(alphaRESULTS[,"P"])~as.factor(alphaRESULTS[,"alpha"]),
  bty="n",xlab=expression(paste(alpha," shape parameter")),
  ylab=expression(paste(Delta," log(L)")),las=1,cex.axis=0.8,
  xlim=c(-0.5,6.5),cex=1.2,pch=21,bg="grey")
abline(h=foo(0.05),lty="dashed",lwd=1,col="red")
boxplot(foo(alphaRESULTS[,"P"])~as.factor(alphaRESULTS[,"alpha"]),
  bty="n",xlab=expression(paste(alpha," shape parameter")),
  ylab=expression(paste(Delta," log(L)")),las=1,cex.axis=0.8,
  xlim=c(-0.5,6.5),add=TRUE,cex=1.2,pch=21,bg="grey")
legend("topright","threshold for statistical significance",
  lwd=1,col="red",cex=0.8,lty="dashed",bty="n")
mtext("B",line=0.5,adj=-0.1,font=2,cex=1.5)

save(alphaRESULTS,rates,trees,x,file="data/alphaRESULTS.rda")

library(vioplot)
vioplot(log(alphaRESULTS[,"est.alpha"])~as.factor(alphaRESULTS[,"alpha"]))
abline(h=log(0.25))
abline(h=log(0.5))
abline(h=log(1))
abline(h=log(2))
abline(h=log(4))  
abline(h=log(8))
