# simulate the data
library("bayesm")
R = 2000
dim = 5
k = 5
sigma = matrix(rep(0.05, dim^2), nrow = dim)
diag(sigma) = 1
sigfac = c(10,20,30,40,50)
mufac = c(1,4,9,16,25)
compsmv = list()
for (i in 1:k){
  compsmv[[i]] = list(mu = mufac[i]*1:dim, sigma = sigfac[i]*sigma)
}

comps = list()
for (i in 1:k) {
  comps[[i]] = list(mu = compsmv[[i]][[1]], rooti = solve(chol(compsmv[[i]][[2]])))
}
pvec = (1:k) / sum(1:k)

nobs = 40
dm = rmixture(nobs, pvec, comps)
Data1 = list(y=dm$x)
ncomp = 9
Prior1 = list(ncomp=ncomp)
out = rnmixGibbs(Data=Data1, Prior=Prior1, Mcmc=Mcmc1)

cat("Summary of Normal Mixture Distribution", fill=TRUE)
summary(out$nmix)

tmom = momMix(matrix(pvec,nrow=1), list(comps))
mat = rbind(tmom$mu, tmom$sd)
cat(" True Mean/Std Dev", fill=TRUE)
print(mat)

## plotting examples
if(1){plot(out$nmix,Data=dm$x)}

