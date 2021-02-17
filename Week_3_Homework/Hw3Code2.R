myGibbs1 = function(x,y,R,keep,BurnIn){
  nvar = ncol(x)   	# the number of independent variables
  nobs = length(y)
  Y = log(y)
  X = c(x[,1],log(x[,2]),x[,3])
  X = matrix(X, nrow = nobs)
  nu = 3
  betabar = c(rep(0, nvar))
  A = diag(rep(0.01, nvar))
  ssq = var(Y)
  XpX = crossprod(X) 	   # t(X)*X
  Xpy = crossprod(X, Y)  # t(X)*y
  sigmasq = as.vector(ssq)
  sigmasqdraw = double(floor((R-BurnIn)/keep)) 
  betadraw = matrix(double(floor( (R-BurnIn) * nvar/keep)), ncol = nvar)
  for (rep in 1:R) {
    
    IR = backsolve(chol(XpX/sigmasq + A), diag(nvar)) 
    # compute the inverse of the root of (XpX/sigmasq+A) 
    btilde = crossprod(t(IR)) %*% (Xpy/sigmasq + A %*% betabar)
    # Different from "runireg, the posterior mean (btilde) 
    # and the posterior variance (crossprod(t(IR))) are conditional on sigmasq

    beta = btilde + IR %*% rnorm(nvar) 
    while (beta[3] <0 ){
      beta = btilde + IR %*% rnorm(nvar)
    }
    res = Y - X %*% beta
    s = t(res) %*% res
    sigmasq = (nu * ssq + s)/rchisq(1, nu + nobs)
    # the scale parameter of the inverted chi-square draw only contain the 
    # the residual of the data. 
    sigmasq = as.vector(sigmasq)
    if (rep > BurnIn){
      if (rep%%keep == 0) {
        mkeep = (rep-BurnIn)/keep
        betadraw[mkeep, ] = beta
        #cat("betaDraw = ", beta, "\n")
        sigmasqdraw[mkeep] = sigmasq
      }
    }
    
  }
  return(list(betadraw = betadraw, sigmasqdraw = sigmasqdraw))
}





myGibbs2 = function(x,y,R,keep,BurnIn){
  nvar = ncol(x)   	# the number of independent variables
  nobs = length(y)
  Y = log(y)
  X = c(x[,1],log(x[,2]),x[,3])
  X = matrix(X, nrow = nobs)
  A = diag(0.01,3)
  nu = 3
  betabar = c(rep(0, nvar))
  sigmatrix = diag(1,3)
  XpX = crossprod(X) 	   # t(X)*X
  Xpy = crossprod(X, Y)  # t(X)*y
  ssq = as.vector(var(Y))
  sigmasq = as.vector(ssq)
  sigmasqdraw = double(floor((R-BurnIn)/keep)) 
  betadraw = matrix(double(floor( (R-BurnIn) * nvar/keep)), ncol = nvar)
  IA = solve(A)
 
  betaone = betabar[c(1,2)]
  betadisplay = betabar[3]
  # cat(V = solve(XpX + A))
  for (rep in 1:R) {
    
    
              
    IR = backsolve(chol(XpX/sigmasq + A), diag(nvar)) 
    V =  (sigmasq * solve(XpX) + solve(A))
    cat("V is: ", V, "\n")
    btilde = crossprod(t(IR)) %*% (Xpy/sigmasq + A %*% betabar)
    betaonetilde = btilde[c(1,2)]
    betadisplaytilde = btilde[3]
    # compute the inverse of the root of (XpX/sigmasq+A) 
    
    # Different from "runireg, the posterior mean (btilde) 
    # and the posterior variance (crossprod(t(IR))) are conditional on sigmasq
    betaonebar = betabar[c(1,2)]
    betadisplaybar = betabar[3]
    V11 = V[1]
    V12 = matrix(V[c(4,7)], nrow = 1)
    V21 = matrix(V[c(2,3)], nrow = 2)
    V22 = matrix(V[c(5,6,8,9)], nrow = 2)
    # cat(betadisplaytilde-solve(rnorm(1,betadisplaytilde-solve(V11)%*%V12%*%(betaone-betaonetilde))), "\n")
    
    cat((solve(V11) %*% V12 %*%(betaonebar - betaonetilde)),betaonebar - betaonetilde ,"\n")

    betadisplay = betadisplaytilde + rtruncnorm(1,-betadisplaytilde,Inf,(betadisplaybar - solve(V11) %*% V12 %*%(betaonebar - betaonetilde)), solve(V11))
    
    cat("betadisplay is: ",betadisplay,"\n")
    
    
    betaone = c(betaonetilde[1] + rnorm(1,(betaonebar-solve(V22)%*%V21%*%(betadisplaybar-betadisplaytilde))[1],solve(V22)[1])
                ,betaonetilde[2]+rnorm(1,(betaonebar-solve(V22)%*%V21%*%(betadisplaybar-betadisplaytilde))[2],solve(V22)[4]) )
    cat("betaone is: ",betaone,"\n")
    beta = matrix(c(betaone,betadisplay),nrow = 3)

    res = Y - X %*% beta
    s = t(res) %*% res
    sigmasq = (nu * ssq + s)/rchisq(1, nu + nobs)
    # the scale parameter of the inverted chi-square draw only contain the 
    # the residual of the data. 
    sigmasq = as.vector(sigmasq)
    # sigmasq = as.vector(3)
    if (rep%%keep == 0 & rep > BurnIn ) {
      mkeep = (rep-BurnIn)/keep
      betadraw[mkeep, ] = beta
      #cat("betaDraw = ", beta, "\n")
      sigmasqdraw[mkeep] = sigmasq
    }
  }
  return(list(betadraw = betadraw, sigmasqdraw = sigmasqdraw))
}



# stimulate the data to test my model
myBetas = c(6.45679,-2.124,0.5535)
obsnum = 1061
myXs = rep(0,obsnum*3)
myYs = rep(0,obsnum)
for (rep in 1:obsnum) {
  myXs[rep] = 1
  myXs[obsnum+rep] = rnorm(1,3,0.25)
  temp = rnorm(1,-0.5,1)
  if (temp > 0){
    myXs[obsnum*2+rep] = temp
  }else{
    myXs[obsnum*2+rep] = 0
  }
  myYs[rep] = exp(myBetas[1] + myBetas[2]*log(myXs[rep+obsnum]) + myBetas[3]*myXs[rep+obsnum*2] + rnorm(1,0,0.1))
}
myXs = matrix(myXs,nrow = obsnum)

result1 = myGibbs1(myXs,myYs,10000,100,5000)
result2 = myGibbs2(myXs,myYs,10000,100,5000)

hist(result1$betadraw[,2])
hist(result2$betadraw[,2])
result = myGibbs2(myXs,myYs,10000,1,5000)
hist(result$betadraw[,3])


data("cheese")
data = cheese[cheese$RETAILER == "ATLANTA - WINN DIXIE",]
Ys = matrix(data$VOLUME)
X = matrix(c(rep(1,61),data[["PRICE"]],data[["DISP"]]),nrow = 61)
result1 = myGibbs1(X,Ys,10000,10,5000)
result2 = myGibbs2(X,Ys,10000,10,5000)

result3 = myGibbs1(X,Ys,100,1,50)
result4 = myGibbs2(X,Ys,100,1,50)
result5 = myGibbs1(X,Ys,100,1,50)
result6 = myGibbs2(X,Ys,100,1,50)
result7 = myGibbs1(X,Ys,100,1,50)
result8 = myGibbs2(X,Ys,100,1,50)

mean(result1$betadraw[,1])
mean(result2$betadraw[,1])
mean(result1$betadraw[,2])
mean(result2$betadraw[,2])
mean(result1$betadraw[,3])
mean(result2$betadraw[,3])
mean(result1$sigmasqdraw)
mean(result2$sigmasqdraw)
hist(result1$sigmasqdraw)
hist(result2$betadraw[,3])
mean(result2$betadraw[,3],)
mean(result4$betadraw[,3],)
mean(result5$betadraw[,3],)
mean(result6$betadraw[,3],)
mean(result7$betadraw[,3],)
mean(result8$betadraw[,3],)