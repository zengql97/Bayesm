lossABeta = c()
for (i in 1: 1000){
  sigmasq = (3*1)/rchisq(1,3)
  A = 0.01*diag(3)
  sqrtInvA = chol(solve(A))
  beta = 0 + as.vector( sqrt(sigmasq) * sqrtInvA %*% rnorm(3))
  lossABeta = c(lossABeta,beta)
}
lossABeta = t (matrix(lossABeta, nrow = 3))
beta1 = lossABeta[,1]
sd(beta1)
beta2 = lossABeta[,2]
sd(beta2)
beta3 = lossABeta[,3]
sd(beta3)


tightABeta = c()

for (i in 1: 1000){
  sigmasq = (3*1)/rchisq(1,3)
  A = 100*diag(3)
  sqrtInvA = chol(solve(A))
  beta = 0 + as.vector( sqrt(sigmasq) * sqrtInvA %*% rnorm(3))
  tightABeta = c(tightABeta,beta)
}
tightABeta = t (matrix(tightABeta, nrow = 3))
beta1 = tightABeta[,1]

sd(beta1)
beta2 = tightABeta[,2]
sd(beta2)
beta3 = tightABeta[,3]
sd(beta3)

boxplot(cbind(lossABeta,tightABeta))