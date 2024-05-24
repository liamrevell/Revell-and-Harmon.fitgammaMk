set.seed(77)

library(phytools)
library(foreach)
library(parallel)
library(doParallel)

niter<-10
nrep<-200
alpha<-0.5
ntaxa<-500

## set number of cores to use
ncores<-min(c(parallel::detectCores()-2,
  niter))
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)

asrRESULTS<-matrix(NA,nrep,8,
  dimnames=list(1:nrep,
    c("delta",
    "gamma.accuracy","est.alpha","gamma.q","gamma.logL",
      "mk.accuracy","mk.q","mk.logL")))

## generate Q matrix for simulation
Q<-matrix(c(
  -1,1,0,
  1,-2,1,
  0,1,-1),3,3,
  dimnames=list(letters[1:3],letters[1:3]))
## design matrix of model
MODEL<-Q
diag(MODEL)<-0

trees<-pbtree(n=ntaxa,nsim=nrep,scale=10)

for(i in 1:nrep){
  ## simulate data
  rates<-rgamma(n=nrow(trees[[i]]$edge),alpha,alpha)
  sim_tree<-trees[[i]]
  sim_tree$edge.length<-trees[[i]]$edge.length*rates
  x<-sim.Mk(sim_tree,Q,internal=TRUE)
  anc<-x[1:sim_tree$Nnode+Ntip(sim_tree)]
  x<-to.matrix(x[sim_tree$tip.label],letters[1:ncol(Q)])
  plotTree(sim_tree,ftype="off",lwd=1,mar=c(1.1,1.1,4.1,1.1))
  tmp<-bquote('N'==.(round(ntaxa))~','~alpha==.(round(alpha,3)))
  title(main=tmp)
  fits<-foreach(j=1:niter)%dopar%{
    result<-NA
    class(result)<-"try-error"
    while(inherits(result,"try-error")){
      result<-try(phytools::fitgammaMk(trees[[i]],x,model=MODEL,
        min.alpha=0.001,pi="fitzjohn",rand_start=TRUE,nrates=8))
    }
    result
  }
  logL<-sapply(fits,logLik)
  gamma_fit<-fits[[which(logL==max(logL))[1]]]
  fits<-foreach(j=1:niter)%dopar%{
    result<-NA
    class(result)<-"try-error"
    while(inherits(result,"try-error")){
      result<-try(phytools::fitMk(trees[[i]],x,model=MODEL,
        pi="fitzjohn",rand_start=TRUE))
    }
    result
  }
  logL<-sapply(fits,logLik)
  mk_fit<-fits[[which(logL==max(logL))[1]]]
  gamma_ace<-ancr(gamma_fit)
  mk_ace<-ancr(mk_fit)
  AA<-to.matrix(anc,letters[1:ncol(Q)])
  mk_accuracy<-1-(1/trees[[i]]$Nnode)*sum(
    (AA-mk_ace$ace)[cbind(1:nrow(AA),apply(AA,1,
      function(x) which(x==1)))])
  gamma_accuracy<-1-(1/trees[[i]]$Nnode)*sum(
    (AA-gamma_ace$ace)[cbind(1:nrow(AA),apply(AA,1,
      function(x) which(x==1)))])
  asrRESULTS[i,]<-c(gamma_accuracy-mk_accuracy,
    gamma_accuracy,gamma_fit$alpha,gamma_fit$rates,
    logLik(gamma_fit),mk_accuracy,mk_fit$rates,logLik(mk_fit))
  if(i>20) ii<-i-19:0
  else ii<-1:i
  print(asrRESULTS[ii,,drop=FALSE])
  cat("\n")
}
parallel::stopCluster(cl=mc)
save(asrRESULTS,file="asrRESULTS.rda")



load("asrRESULTS.rda")
prop_diff<-asrRESULTS[,"delta"]/asrRESULTS[,"mk.accuracy"]
hist(prop_diff,breaks=seq(-0.15,0.15,by=0.01),main="",
  xlab=expression(paste(Delta," reconstruction accuracy")),
  ylab="frequency",las=1,cex.axis=0.8,ylim=c(0,40))
grid()
clip(par()$usr[1],par()$usr[2],0,par()$usr[4])
abline(v=0,lwd=2)
tmp<-bquote(Gamma ~ .("model more accurate") %->% "")
text(x=0,y=40,tmp,pos=4,cex=0.8)
tmp<-bquote("" %<-% "Mk model more accurate")
text(x=0,y=40,tmp,pos=2,cex=0.8)
