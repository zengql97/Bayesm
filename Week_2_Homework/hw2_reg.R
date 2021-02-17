function(y, X, betabar, A, nu, ssq){
  n = length(y)
  # n is the number of observasions
  k = ncol(X)
  # k is the number of independent variables. (k-1,X[1] is 1, for beta_0)
  RA = chol(A)
  # RA will be used later
  W = rbind(X,RA)
  # W = (X; RA)
  z = c(y, as.vector(RA %*% betabar))
  # z = (y, ra*betabar)
  IR = backsolve(chol(crossprod(W)), diag(k))
  # croosprod(W) = X'X + A
  # IR is the normalized inverse of X'X + A
  # IR * IR' = (X'X+A)^-1
  btilde = crossprod(t(IR)) %*% crossprod(W,z)
  # t(W)* z = (x'y + A*betahat)
  # btilde = (X'X+A)-1 * (x'y + A*betahat)
  res = z - W %*% btilde
  # W x btilde = (yhat,ra*btilde)
  s = crossprod(res)
  # s = vs^2
  sigmasq = (nu*ssq + s)/rchisq(1,nu+n)
  # draw sig = v1s1^2/chi-sq(v_1)
  beta = btilde + as.vector(sqrt(sigmasq)) * IR %*% rnorm(k)
  # draw beta = N(betatlide,sig(X'X+A)^-1)
  list(beta = beta,sigmasq = sigmasq)
  
}