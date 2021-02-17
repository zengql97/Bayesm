library(MASS, lib.loc = "C:/Program Files/R/R-4.0.3/library")
setwd("D:\\Dropbox\\SMU\\21Spring\\MKTG707\\Week_4_Homework")
binaryProbitData = function(delta,Obs = 1000 ,Nx = 1,meanX = 0,sig = 1){
  if (Nx == 1){
    X = rnorm(Obs,meanX,sigma)
  }else{
    X = mvrnorm( n =Obs,mu = meanX, Sigma = diag(sig,3) )
  }
  
  Z = X%*%delta + rnorm(Obs,0,100)
  Y = as.integer(Z > 0)
  return(list (X = X, y = Y) )
}
Data  = binaryProbitData(c(-1,-2,3),Obs = 1000, Nx = 3, meanX = c(1,2,3), sig = c(1,2,3) )
binaryProbitData = function(delta,Obs = 1000 ,Nx = 1,meanX = 0,sig = 1){
  if (Nx == 1){
    X = rnorm(Obs,meanX,sigma)
  }else{
    X = mvrnorm( n =Obs,mu = meanX, Sigma = diag(sig,3) )
  }
  
  Z = X%*%delta + rnorm(Obs,0,100)
  Y = as.integer(Z > 0)
  return(list (X = X, y = Y) )
}





simmnl= function(X,p,n,beta) {
  #   note: create X array with 2 alt.spec vars
  k=length(beta)

  Xbeta=X%*%beta # now do probs
  p=nrow(Xbeta)/n
  Xbeta=matrix(Xbeta,byrow=TRUE,ncol=p)
  Prob=exp(Xbeta)
  iota=c(rep(1,p))
  denom=Prob%*%iota
  Prob=Prob/as.vector(denom)
  # draw y
  y=vector("double",n)
  ind=1:p
  for (i in 1:n) 
  { yvec=rmultinom(1,1,Prob[i,]); y[i]=ind%*%yvec }
  return(list(y=y,X=X,beta=beta,prob=Prob))
}

simmnp = function(X,p, n, beta, sigma) {
  k=length(beta)
  indmax = function(x) {which(max(x)==x)}

  Xbeta = X%*%beta
  w = as.vector(crossprod(chol(sigma),matrix(rnorm((p-1)*n),ncol=n))) + Xbeta
  w = matrix(w, ncol=(p-1), byrow=TRUE)
  maxw = apply(w, 1, max)
  y = apply(w, 1, indmax)
  y = ifelse(maxw < 0, p, y)
  return(list(y=y, X=X, beta=beta, sigma=sigma))
}

