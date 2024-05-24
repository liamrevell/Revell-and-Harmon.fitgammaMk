set.seed(88)

library(foreach)
library(parallel)
library(doParallel)

lr.test<-function(obj1,obj2){
  l1<-logLik(obj1)
  l2<-logLik(obj2)
  LR<-2*(l2-l1)
  as.numeric(
    pchisq(LR,df=attr(l2,"df")-attr(l1,"df"),
      lower.tail=FALSE))
}

niter<-10
nrep<-100
alpha<-0.5
ntaxa<-c(50,100,200,400,800,1600)

## set number of cores to use
ncores<-min(c(parallel::detectCores()-2,
  niter))
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

powerRESULTS<-matrix(NA,nrep*length(ntaxa),5,
  dimnames=list(1:(nrep*length(ntaxa)),
    c("N","est.alpha","est.q","logL","P")))
## generate Q matrix for simulation
Q<-matrix(c(
  -1,1,0,
  1,-2,1,
  0,1,-1),3,3,
  dimnames=list(letters[1:3],letters[1:3]))
## design matrix of model
MODEL<-Q
diag(MODEL)<-0
for(i in 1:length(ntaxa)){
  ## generate random tree(s)
  trees<-pbtree(n=ntaxa[i],scale=1,nsim=nrep)
  for(j in 1:nrep){
    ## simulate data
    rates<-rgamma(n=nrow(trees[[j]]$edge),alpha,alpha)
    sim_tree<-trees[[j]]
    sim_tree$edge.length<-trees[[j]]$edge.length*rates
    x<-sim.Mk(sim_tree,Q)
    x<-to.matrix(x,letters[1:ncol(Q)])
    plotTree(sim_tree,ftype="off",lwd=1,mar=c(1.1,1.1,4.1,1.1))
    tmp<-bquote('N'==.(round(ntaxa[i]))~','~alpha==.(round(alpha,3)))
    title(main=tmp)
    fits<-foreach(k=1:niter)%dopar%{
      result<-NA
      class(result)<-"try-error"
      while(inherits(result,"try-error")){
        result<-try(phytools::fitgammaMk(trees[[j]],x,model=MODEL,
          min.alpha=0.001,pi="fitzjohn",rand_start=TRUE,nrates=8))
      }
      result
    }
    logL<-sapply(fits,logLik)
    best_fit<-fits[[which(logL==max(logL))[1]]]
    fits<-foreach(k=1:niter)%dopar%{
      result<-NA
      class(result)<-"try-error"
      while(inherits(result,"try-error")){
        result<-try(phytools::fitMk(trees[[j]],x,model=MODEL,
          pi="fitzjohn",rand_start=TRUE))
      }
      result
    }
    logL<-sapply(fits,logLik)
    fit_null<-fits[[which(logL==max(logL))[1]]]
    powerRESULTS[(i-1)*nrep+j,]<-
      c(ntaxa[i],best_fit$alpha,best_fit$rates,
        logLik(best_fit),lr.test(fit_null,best_fit))
    if((((i-1)*nrep+j))>20) ii<-((i-1)*nrep+j)-19:0
    else ii<-1:((i-1)*nrep+j)
    print(round(powerRESULTS[ii,],5))
    cat("\n")
  }
}
stopCluster(cl=mc)
save(powerRESULTS,file="powerRESULTS.rda")

load("powerRESULTS.rda")
ntaxa<-c(50,100,200,400,800,1600)
alpha<-0.5
par(bty="n",mfrow=c(1,2),mar=c(5.1,4.1,2.1,1.1))
boxplot(powerRESULTS[,"est.alpha"]~as.factor(powerRESULTS[,"N"]),log="y",
  bty="n",xlab="number of taxa",
  ylab=expression(paste("estimated value of ",alpha)),las=1,cex.axis=0.8,
  xlim=c(-0.5,6.5),ylim=c(0.001,1000),axes=FALSE)
abline(h=alpha,lty="dashed",col="blue")
boxplot(powerRESULTS[,"est.alpha"]~as.factor(powerRESULTS[,"N"]),log="y",
  add=TRUE,cex=1.2,pch=21,bg="grey",las=1,cex.axis=0.8,axes=FALSE)
axis(1,at=1:6,labels=ntaxa,las=1,cex.axis=0.8)
axis(2,at=10^(-3:3),labels=c(0.001,0.01,0.1,1.0,10.0,100,1000),
  las=1,cex.axis=0.8)
tmp<-bquote(alpha==.(round(alpha,3)))
text(x=-0.7,y=1.4*alpha,tmp,pos=4,cex=0.8)
legend("bottomright",expression(paste("generating value of ",alpha)),
  lwd=1,col="blue",lty="dashed",bty="n",cex=0.8)
mtext("A",line=0.5,adj=-0.1,font=2,cex=1.5)
PP<-sapply(ntaxa,function(n,Y) Y[which(Y[,"N"]==n),"P"],
  Y=powerRESULTS)
power<-apply(PP,2,function(x) mean(x<=0.05))
plot(1:6,power,type="b",pch=21,cex=1.2,bg="grey",
  xlim=c(-0.5,6.5),ylim=c(0,1),las=1,cex.axis=0.8,
  axes=FALSE,xlab="number of taxa",ylab="estimated power")
axis(2,las=1,cex.axis=0.8)
axis(1,at=1:6,labels=ntaxa,las=1,cex.axis=0.8)
abline(h=0.05,lty="dashed",col="blue")
text(x=-0.7,y=0.08,"P = 0.05",pos=4,cex=0.8)
mtext("B",line=0.5,adj=-0.1,font=2,cex=1.5)

