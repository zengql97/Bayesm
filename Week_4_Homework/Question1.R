Data  = binaryProbitData(c(-1,-2,3),Obs = 1000, Nx = 3, meanX = c(1,2,3), sig = c(1,2,3) )
Mcmc = list(R = 10000,keep = 1)
Result1 = myrbprobitGibbs(Data = Data, Mcmc = Mcmc)

Result2 = rbprobitGibbs(Data = Data, Mcmc = Mcmc)

 


mean(Result2$betadraw[,1])
mean(Result1$betadraw[,1]/(Result1$sigmadraw^0.5))
plot(Result1$betadraw[,1])
