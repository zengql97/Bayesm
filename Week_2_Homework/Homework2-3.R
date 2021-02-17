function(Y, X, Bbar, A, nu, V){
  # n is the number of observations
  n = nrow(Y)
  # m is the number of dependent variables
  m = ncol(Y)
  # k is the number of 
  k = ncol(X)
  # RA is the chol decomposition of A
  RA = chol(A)
  # W = (X; RA)
  W = rbind(X, RA)
  # Z = (Y; RA %*% Bbar)
  Z = rbind(Y, RA%*% Bbar)
  # IR is the normalized inverse matrix W
  IR = backsolve(chol(crossprod(W)), diag(k))
  # The Max likelihood B for the posterior
  Btilde = crossprod(t(IR)) %*% crossprod(W,Z)
  # IRIR'(W'Z) = (X'X+A)^-1(X'Y + ABbar)
  # S is the residual between y, y_hat and RA Bbar, RA*btilde
  S = crossprod(Z-W %*% Btilde)
  # draw Sigma from Inverted Wishart Distribution, nu = nu + n, 
  # V+S is inverted to draw the IW Prior.
  # n <- degree of freedom
  rwout = rwishart(nu+n, chol2inv(chol(V+S)))
  B = Btilde + IR %*% matrix(rnorm(m*k), ncol = m) %*% t(rwout$CI)
  return(list(B = B, Sigma = rwout$IW ))
  #B = draw of regression coefficient matrix
  
  
  
  
}




if(nchar(Sys.getenv("LONG_TEST")) != 0) {R=2000} else {R=10}
set.seed(66)
R = 20000
n =200
m = 2
X = cbind(rep(1,n),runif(n))
k = ncol(X)
B = matrix(c(1,2,-1,3), ncol=m)
Sigma = matrix(c(1, 0.5, 0.5, 1), ncol=m)
RSigma = chol(Sigma)
Y = X%*%B + matrix(rnorm(m*n),ncol=m)%*%RSigma

betabar = rep(0,k*m)
Bbar = matrix(betabar, ncol=m)
A = diag(rep(0.01,k))
nu = 3
V = nu*diag(m)

betadraw = matrix(double(R*k*m), ncol=k*m)
Sigmadraw = matrix(double(R*m*m), ncol=m*m)

for (rep in 1:R) {
  out = rmultireg(Y, X, Bbar, A, nu, V)
  betadraw[rep,] = out$B
  Sigmadraw[rep,] = out$Sigma
}

cat(" Betadraws ", fill=TRUE)
mat = apply(betadraw, 2, quantile, probs=c(0.01, 0.05, 0.5, 0.95, 0.99))
mat = rbind(as.vector(B),mat)
rownames(mat)[1] = "beta"
print(mat)

cat(" Sigma draws", fill=TRUE)
mat = apply(Sigmadraw, 2 ,quantile, probs=c(0.01, 0.05, 0.5, 0.95, 0.99))
mat = rbind(as.vector(Sigma),mat); rownames(mat)[1]="Sigma"
print(mat)