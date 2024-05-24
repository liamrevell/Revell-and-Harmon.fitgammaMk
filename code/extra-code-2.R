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
niter<-10
nrep<-20
alpha<-c(0.25,0.5,1,2,4,8)
ntaxa<-2000
## set number of cores to use
ncores<-min(c(parallel::detectCores()-2,
  niter))
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)
## results matrix
alphaRESULTS<-matrix(NA,nrep*length(alpha),5,
  dimnames=list(1:(nrep*length(alpha)),
    c("alpha","est.alpha","est.q","logL","P")))
## generate random tree(s)
trees<-pbtree(n=ntaxa,scale=1,nsim=nrep)
## generate Q matrix for simulation
Q<-matrix(c(
  -1,1,0,
  1,-2,1,
  0,1,-1),3,3,
  dimnames=list(letters[1:3],letters[1:3]))
## design matrix of model
MODEL<-Q
diag(MODEL)<-0
for(i in 1:length(alpha)){
  for(j in 1:nrep){
    ## simulate data
    rates<-rgamma(n=nrow(trees[[j]]$edge),alpha[i],alpha[i])
    sim_tree<-trees[[j]]
    sim_tree$edge.length<-trees[[j]]$edge.length*rates
    x<-sim.Mk(sim_tree,Q)
    x<-to.matrix(x,letters[1:ncol(Q)])
    plotTree(sim_tree,ftype="off",lwd=1,mar=c(1.1,1.1,4.1,1.1))
    tmp<-bquote(alpha==.(round(alpha[i],3)))
    title(main=tmp)
    fits<-foreach(k=1:niter)%dopar%{
      phytools::fitgammaMk(trees[[j]],x,model=MODEL,min.alpha=0.001,
        pi="fitzjohn",rand_start=TRUE,nrates=8)
    }
    logL<-sapply(fits,logLik)
    best_fit<-fits[[which(logL==max(logL))[1]]]
    fits<-foreach(k=1:niter)%dopar%{
      phytools::fitMk(trees[[j]],x,model=MODEL,
        pi="fitzjohn",rand_start=TRUE)
    }
    logL<-sapply(fits,logLik)
    fit_null<-fits[[which(logL==max(logL))[1]]]
    alphaRESULTS[(i-1)*nrep+j,]<-
      c(alpha[i],best_fit$alpha,best_fit$rates,
        logLik(best_fit),lr.test(fit_null,best_fit))
    ## print(alphaRESULTS[1:((i-1)*nrep+j),])
  }
}
stopCluster(cl=mc)
save(alphaRESULTS,file="alphaRESULTS.rda")