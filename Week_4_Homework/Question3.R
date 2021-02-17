fsh=function() 
{
  # 
  # P. Rossi
  # revision history: 3/27/05
  #
  # Purpose:
  #  function to flush console (needed only under windows)
  #
  if (Sys.info()[1] == "Windows") flush.console()
  return()
}
n=200; p=4; beta=c(10,10,10,15,5,9)
simmnl= function(p,n,beta) {
  #   note: create X array with 2 alt.spec vars
  k=length(beta)
  X1=matrix(runif(n*p,min=-1,max=1),ncol=p)
  X2=matrix(runif(n*p,min=-1,max=1),ncol=p)
  X3=matrix(runif(n*p,min=-1,max=1),ncol=p)
  X=createX(p,na=k-p+1,nd=NULL,Xd=NULL,Xa=cbind(X1,X2,X3),base=1)
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
simout=simmnl(p,n,beta)
Data = list(X = simout$X, y = simout$y, p = 4)
Mcmc = list(R = 100000, keep = 1)
r1 = rmnlIndepMetrop(Data = Data, Mcmc = Mcmc)
r2 = rmnlRWMetrop(Data = Data, Mcmc = Mcmc)
plot(r2$betadraw[,2])


n=200; p=3; beta=c(0.14,0.4,1.5,1.6,1.7)
simout=simmnp(p,n,beta,sigma = diag(1,2))
Data = list(X = simout$X, y = simout$y, p = 3)
Mcmc = list(R = 10000, keep = 10)
r2 = rmnpGibbs(Data = Data, Mcmc = Mcmc)
plot(r2$betadraw)
n=200; p=3; beta=c(0.14,0.4,1.5,1.6,1.7)
simout=simmnl(p,n,beta)
Data = list(X = simout$X, y = simout$y, p = 3)
Mcmc = list(R = 10000, keep = 10)

r1 = rmnlIndepMetrop(Data = Data, Mcmc = Mcmc)

plot(r2$betadraw[,4])



beta=c(0.14,0.4,1.5,1.6,1.7)

k=length(beta)
X1=matrix(runif(n*p,min=-1,max=1),ncol=p)
X2=matrix(runif(n*p,min=-1,max=1),ncol=p)
X3=matrix(runif(n*p,min=-1,max=1),ncol=p)
XA=createX(p,na=k-p+1,nd=NULL,Xd=NULL,Xa=cbind(X1,X2,X3),base=p)
XB=createX(p,na=k-p+1,nd=NULL,Xd=NULL,Xa=cbind(X1,X2,X3),DIFF=TRUE,base=p)

simout1=simmnp(XB,p,n,beta,sigma = diag(1,2))
simout2=simmnl(XA,p,n,beta)
Mcmc = list(R = 10000, keep = 10)
Data1 = list(X = simout1$X, y = simout1$y, p = 3)
Data2 = list(X = simout2$X, y = simout2$y, p = 3)
r1 = rmnpGibbs(Data = Data1, Mcmc = Mcmc)
r2 = rmnlIndepMetrop(Data = Data2, Mcmc = Mcmc)
plot(r4$betadraw[,4])


Data3 = list(X = simout1$X, y = simout2$y, p = 3)
r3 = rmnpGibbs(Data = Data3, Mcmc = Mcmc)
Data4 = list(X = simout2$X, y = simout1$y, p = 3)
r4 = rmnlIndepMetrop(Data = Data4, Mcmc = Mcmc)