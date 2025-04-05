## load packages
library(phytools)
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
nrep<-100
alpha<-c(0.25,2)
ntaxa<-c(50,100,200,400,800,1600)
## create data matrix for results
powerRESULTS<-matrix(NA,nrep*length(ntaxa)*length(alpha),7,
  dimnames=list(1:(nrep*length(ntaxa)*length(alpha)),
    c("alpha","N","est.alpha","est.q","logL","logL(null)","P")))
## set number of cores to use
ncores<-min(c(parallel::detectCores()-2,
  niter))
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)
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
for(i in 1:length(alpha)){
  for(j in 1:length(ntaxa)){
    ## generate random tree(s)
    set.seed(j) ## set seed internally (this is on purpose)
    trees<-pbtree(n=ntaxa[j],scale=10,nsim=nrep)
    for(k in 1:nrep){
      ## simulate data
      rates<-rgamma(n=nrow(trees[[k]]$edge),alpha[i],alpha[i])
      sim_tree<-trees[[k]]
      sim_tree$edge.length<-trees[[k]]$edge.length*rates
      x<-sim.Mk(sim_tree,Q)
      x<-to.matrix(x,letters[1:ncol(Q)])
      plotTree(sim_tree,ftype="off",lwd=1,mar=c(1.1,1.1,4.1,1.1))
      tmp<-bquote('N'==.(round(ntaxa[j]))~','~alpha==.(round(alpha[i],3)))
      title(main=tmp)
      fits<-foreach(l=1:niter)%dopar%{
        result<-NA
        class(result)<-"try-error"
        while(inherits(result,"try-error")){
          result<-try(
            R.utils::withTimeout({
              phytools::fitgammaMk(trees[[k]],x,model=MODEL,
                min.alpha=0.001,pi="fitzjohn",rand_start=TRUE,nrates=8)
            },timeout=1200,onTimeout="error"),silent=TRUE)
        }
        result
      }
      logL<-sapply(fits,logLik)
      best_fit<-fits[[which(logL==max(logL))[1]]]
      fits<-foreach(l=1:niter)%dopar%{
        result<-NA
        class(result)<-"try-error"
        while(inherits(result,"try-error")){
          result<-try(
            R.utils::withTimeout({
              phytools::fitMk(trees[[k]],x,model=MODEL,
                pi="fitzjohn",rand_start=TRUE)
            },timeout=1200,onTimeout="error"),silent=TRUE)
        }
        result
      }
      logL<-sapply(fits,logLik)
      fit_null<-fits[[which(logL==max(logL))[1]]]
      powerRESULTS[(i-1)*(nrep*length(ntaxa))+(j-1)*nrep+k,]<-
        c(alpha[i],ntaxa[j],best_fit$alpha,best_fit$rates,
          logLik(best_fit),logLik(fit_null),
          lr.test(fit_null,best_fit))
        if(((i-1)*(nrep*length(ntaxa))+(j-1)*nrep+k)>20) 
          ii<-((i-1)*(nrep*length(ntaxa))+(j-1)*nrep+k)-19:0
        else ii<-1:((i-1)*(nrep*length(ntaxa))+(j-1)*nrep+k)
      cat("\nCurrent time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
      print(round(powerRESULTS[ii,],5))
      cat("\n")
    }
  }
}
stopCluster(cl=mc)
save(powerRESULTS,file="data/powerRESULTS.rda")
