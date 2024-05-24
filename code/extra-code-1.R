ncores<-10
## create parallel socket cluster
mc<-makeCluster(ncores,type="PSOCK")
registerDoParallel(cl=mc)
a1<-a2<-vector()
alpha<-0.25
for(i in 1:10){
  tree<-pbtree(n=800,scale=10)
  q<-0.2
  Q<-matrix(c(
    -q,q,0,
    q,-2*q,q,
    0,q,-q),3,3,
    dimnames=list(letters[1:3],letters[1:3]))
  rates<-rgamma(n=nrow(tree$edge),alpha,alpha)
  sim_tree<-tree
  sim_tree$edge.length<-tree$edge.length*rates
  cat("Edge lengths in original & rate heterogenous tree:\n")
  print(sum(tree$edge.length))
  print(sum(sim_tree$edge.length))
  cat("------------\n\n")
  x1<-sim.Mk(tree,Q)
  x2<-sim.Mk(sim_tree,Q)
  dev.off()
  h<-max(nodeHeights(tree))
  plotTree(tree,ftype="off",lwd=1,
    xlim=c(0,1.1)*h,col="grey")
  cols<-setNames(palette()[1:3],levels(x))
  for(j in 1:length(x)){
    polygon(c(1,1.05,1.05,1)*h,j+c(-0.5,-0.5,0.5,0.5),col=cols[x1[j]],
      border=FALSE)
    polygon(c(1.05,1.1,1.1,1.05)*h,j+c(-0.5,-0.5,0.5,0.5),col=cols[x2[j]],
      border=FALSE)
  }
  fits<-foreach(k=1:10)%dopar%{
    phytools::fitMk(tree,x1,model=MODEL,pi="fitzjohn",rand_start=TRUE)
  }
  fit_mk1<-fits[[which(sapply(fits,logLik)==max(sapply(fits,logLik)))[1]]]
  fits<-foreach(k=1:10)%dopar%{
    phytools::fitgammaMk(tree,x1,model=MODEL,min.alpha=0.001,
      pi="fitzjohn",rand_start=TRUE,nrates=8)
  }
  fit_gamma1<-fits[[which(sapply(fits,logLik)==max(sapply(fits,logLik)))[1]]]
  cat("Model comparison for the Mk data:\n")
  anova(fit_mk1,fit_gamma1)
  a1[i]<-fit_gamma1$alpha
  print(a1[i])
  cat("--------------\n\n")
  cat("Model comparison for the Mk+G data:\n")
  fit_mk2<-fitMk(tree,x2,model=MODEL,pi="fitzjohn",rand_start=TRUE,)
  fit_gamma2<-fitgammaMk(tree,x2,model=MODEL,min.alpha=0.001,
    pi="fitzjohn",rand_start=TRUE,nrates=8)
  anova(fit_mk2,fit_gamma2)
  a2[i]<-fit_gamma2$alpha
  print(a2)
  cat("--------------\n\n")
}
stopCluster(cl=mc)
