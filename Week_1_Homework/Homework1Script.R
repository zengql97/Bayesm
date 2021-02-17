library(dplyr)
library(bayesm)
Rep <- 100
beta = c(1,2,3,4,5)
Xs = matrix(double(Rep*5),ncol = 5)
Draws = matrix(double(Rep*1),ncol = 1)
for (rep in 1:Rep) {Xs[rep,]=c(1,runif(4,0,5))}
for (rep in 1:Rep) {Draws[rep,]=sum(beta*Xs[rep,])+rnorm(1,0,1)}
# ggplot() + geom_point(mapping = aes(x= Xs[,1], y = Draws) )+theme_light() +geom_smooth(aes(x = Xs[,1], y = Draws),method = 'glm') + xlab("X[1]")
# ggplot() + geom_point(mapping = aes(x= Xs[,2], y = Draws) )+theme_light() +geom_smooth(aes(x = Xs[,2], y = Draws),method = 'glm') + xlab("X[2]")
# ggplot() + geom_point(mapping = aes(x= Xs[,3], y = Draws) )+theme_light() +geom_smooth(aes(x = Xs[,3], y = Draws),method = 'glm') + xlab("X[3]")
# ggplot() + geom_point(mapping = aes(x= Xs[,4], y = Draws) )+theme_light() +geom_smooth(aes(x = Xs[,4], y = Draws),method = 'glm') + xlab("X[4]")
y = Draws
X = Xs
myreg = function(y,X) {
  Rep = dim(X)[1]
  Length = dim(X)[2]

  bhat = solve(t(X)%*%X[,]) %*% t(X) %*%y
  yhat = matrix(double(Rep*1),ncol = 1)
  for (rep in 1:Rep) {yhat[rep,]=sum(bhat*Xs[rep,])}
  resid = yhat - y
  ybar = mean(y)
    SST = sum( (y - ybar)^2 )
    SSE = sum( (resid)^2 )
    SSR = sum( (yhat - ybar)^2 )
    Rsq = SSR/SST
    sigma = SSE/(Rep - Length - 1)
    stderror = sqrt(SSE/(Rep - Length - 1))
    covb =sigma * solve(t(X)%*%X[,])
    seb = matrix(double(Length*1),ncol = 1)
    for (i in 1:Length) {seb[i] = sqrt(covb[i,i])}
    tstat = bhat/seb
    list(bhat=bhat, yhat=yhat, resid=resid, ybar=ybar, SST=SST, SSE=SSE, SSR=SSR, Rsq=Rsq, stderror=stderror, seb=seb, tstat=tstat,sr = sum(resid))
}
myreg(y,X)
data("cheese")
Retailer <- as.vector(cheese[[1]])
Retailer[1:10].unique()
filter(cheese, RETAILER %in% Retailer[1:10])
SelectedCheese = filter(cheese, RETAILER %in% Retailer[1:10])
Y <- log(as.vector(SelectedCheese[[2]]))
Xs <- c(rep(1,631),as.vector(SelectedCheese[[3]]), as.vector(SelectedCheese[[4]]))
Xs <- matrix(Xs, ncol = 3)

myreg(Y,Xs)
