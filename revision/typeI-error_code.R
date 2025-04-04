library(phytools)

set.seed(99)
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
niter<-2 #10
nrep<-2 #100
ntaxa<-c(50,100,200,400,800,1600)
## set number of cores to use
ncores<-min(c(parallel::detectCores()-2,
  niter))
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)
typeIerrorRESULTS<-matrix(NA,
  nrep*(length(ntaxa)+1),6,
  dimnames=list(1:(nrep*(length(ntaxa)+1)),
    c("N","est.alpha","est.q","logL","logL(null)","P")))
## generate Q matrix for simulation
Q<-matrix(c(
  -0.5, 0.5, 0.0,
  0.5,-1.0, 0.5,
  0.0, 0.5,-0.5),3,3,
  dimnames=list(letters[1:3],letters[1:3]))
## design matrix of model
MODEL<-matrix(c(
  0,1,0,
  1,0,1,
  0,1,0),3,3,
  dimnames=list(letters[1:3],letters[1:3]))
for(i in 1:length(ntaxa)){
  ## generate random tree(s)
  set.seed(i) ## set seed internally (this is on purpose)
  trees<-pbtree(n=ntaxa[i],scale=10,nsim=nrep)
  for(j in 1:nrep){
    ## simulate data
    x<-sim.Mk(trees[[j]],Q)
    x<-to.matrix(x,letters[1:ncol(Q)])
    plotTree(trees[[j]],ftype="off",lwd=1,mar=c(1.1,1.1,4.1,1.1))
    tmp<-bquote('N'==.(round(ntaxa[j]))*","~alpha==.(quote(infinity)))
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
    typeIerrorRESULTS[(i-1)*nrep+j,]<-
      c(ntaxa[i],best_fit$alpha,best_fit$rates,
          logLik(best_fit),logLik(fit_null),
          lr.test(fit_null,best_fit))
    if(((i-1)*nrep+j)>20) ii<-((i-1)*nrep+j)-19:0
    else ii<-1:((i-1)*nrep+j)
    print(round(typeIerrorRESULTS[ii,,drop=FALSE],5))
    cat("\n")
  }
}
stopCluster(cl=mc)
## save(typeIerrorRESULTS,file="data/typeIerrorRESULTS.rda")
