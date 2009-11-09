### Internal functions for
### mv.2way.est
###



###
###
###

`estimate.identity` <- function(x, block, treatment)
{
  d <- dim(x)[2]
  k <- nlevels(treatment)
  n <- nlevels(block)
  N <- dim(x)[1]

  #Calculate the treatment difference estimates
  res<-vector("list", k-1)
  for(i in 2:k){
    diff<-x[treatment==levels(treatment)[i],]-
      x[treatment==levels(treatment)[1],]
    delta<-apply(diff,2,mean)
    METHOD<-"affine equivariant treatment difference estimate (identity score)"
    dn<-paste("Delta",levels(treatment)[i],"-Delta",levels(treatment)[1],sep="")
    S<-(1/N)*t(x)%*%x
    res[[i-1]]<-list(location=delta,vcov=S,est.name=METHOD,dname=dn)
    class(res[[i-1]])<-"mvloc"
  }

  return(res)
}


###
###
###

`estimate.sign` <- function(x, block, treatment,
                          stand=c("outer","inner"),
                          eps=1.0e-15,maxiter=100)
{
  d <- dim(x)[2]
  k <- nlevels(treatment)
  n <- nlevels(block)
  N <- dim(x)[1]

  switch(stand,
         "inner"=
         {
           METHOD <- "affine equivariant treatment difference estimate (sign score)"
           sqrtG<-inner.sign(x=x,block=block,eps=eps,maxiter=maxiter)$sqrtG
           z <- x%*%t(solve(sqrtG))    
           Delta<-vector("list", k-1)
           for(i in 2:k){
              delta<-spatial.median(z[treatment==levels(treatment)[i],],eps=eps,maxiter=maxiter)-
                     spatial.median(z[treatment==levels(treatment)[1],],eps=eps,maxiter=maxiter)
              Delta[[i-1]]<-delta%*%t(sqrtG)
           }
         },
         "outer"=
         {
           METHOD <- "treatment difference estimate (sign score)"
           Delta<-vector("list", k-1)
           for(i in 2:k){
              delta<-spatial.median(x[treatment==levels(treatment)[i],],eps=eps,maxiter=maxiter)-
                     spatial.median(x[treatment==levels(treatment)[1],],eps=eps,maxiter=maxiter)
              Delta[[i-1]]<-delta
           } 
         }
         )

  res<-vector("list", k-1) 
  S<-NULL
  for(i in 2:k){
    dn<-paste("Delta",levels(treatment)[i],"-Delta",levels(treatment)[1],sep="")
    res[[i-1]]<-list(location=Delta[[i-1]],vcov=S,est.name=METHOD,dname=dn)
    class(res[[i-1]])<-"mvloc"
  }
  
  return(res)
}



###
###
###

`estimate.rank` <- function(x, block, treatment,
                            stand=c("outer","inner"),
                            eps=1.0e-15,maxiter=100)
{
  d <- dim(x)[2]
  k <- nlevels(treatment)
  n <- nlevels(block)
  N <- dim(x)[1]        

  Delta<-vector("list", k)
  for(ii in 1:k){
     Delta[[ii]]<-vector("list", k)
     Delta[[ii]][[ii]]<-matrix(0,1,d)
  }

  #Calculate the unadjusted treatment difference estimates
  switch(stand,
         "inner"=
         {
           METHOD <- "affine equivariant treatment difference estimate (rank score)"
           rnk<-inner.rank(x=x,block=block,eps=eps,maxiter=maxiter)
           rx<-rnk$rx
           sqrtG<-rnk$sqrtG
           z <- x%*%t(solve(sqrtG))    
           for(i in 1:(k-1)){
             for(j in (i+1):k){
                delta<-spatial.median(z[treatment==levels(treatment)[j],]-
                                      z[treatment==levels(treatment)[i],],eps=eps,maxiter=maxiter)
                Delta[[i]][[j]]<-delta%*%t(sqrtG)
                Delta[[j]][[i]]<- -Delta[[i]][[j]]
             }
           }
         },
         "outer"=
         {
           METHOD <- "treatment difference estimate (rank score)"
           rx<-outer.rank(x=x,block=block)
           for(i in 1:(k-1)){
             for(j in (i+1):k){
              delta<-spatial.median(x[treatment==levels(treatment)[j],]-
                                    x[treatment==levels(treatment)[i],],eps=eps,maxiter=maxiter)
              Delta[[i]][[j]]<-delta
              Delta[[j]][[i]]<- -Delta[[i]][[j]]
             }
           } 
         }
         )

  #Calculate the asymptotic covariance matrix of the location estimates
  B<-(1/N)*t(rx)%*%rx
  A<-matrix(0,d,d)
  for(i in 1:(k-1)){
    for(j in (i+1):k){
      diff<- x[treatment==levels(treatment)[j],]-
             x[treatment==levels(treatment)[i],]
      r<-SpatialNP:::norm(diff)
      u<-spatial.sign(diff,center=FALSE,shape=FALSE)
      A<-A+(1/n)*(sum(1/r)*diag(d)-t(u)%*%diag(1/r)%*%u)
    }
  }
  A<-(2/(k*(k-1)))*A
  inv.A<-solve(A)
  S<-(2*k/(k-1))*inv.A%*%B%*%inv.A


  #Calculate the adjusted treatment difference estimates
  res<-vector("list", k-1)
  for(i in 2:k){
    dn<-paste("Delta",levels(treatment)[i],"-Delta",levels(treatment)[1],sep="")
    s1<-rep(0,d)
    s2<-rep(0,d)
    for(ii in 1:k){
      s1<-s1+Delta[[1]][[ii]]
      s2<-s2+Delta[[i]][[ii]]
    }
    delta<-(s1-s2)/k
    res[[i-1]]<-list(location=delta,vcov=S,est.name=METHOD,dname=dn)
    class(res[[i-1]])<-"mvloc"
  }
    
  return(res)
}
