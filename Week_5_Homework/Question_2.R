# import the data from "bayesm"

library("bayesm")
data("bank")
# define the function to map -1 to 2
ZeroToTwo = function(x){
  result = as.integer(x)
  if (x == 0){
    result = 2
  }
  return(result)
}
# define X of the estimation
X = as.data.frame(bank$choiceAtt)
# the set of users
IDlist = unique(X$id)
# prepare to split the data x, based on the ID
lgtdata = list()
# init lgt data

for (i in IDlist){
  suby = as.integer(map(X[X$id == i,][[2]],ZeroToTwo))
  # the dependent var
  subX = as.matrix(X[ X$id == i,][c(3:16)])
  # independent vat
  nobs = length(subX)/14
  tempX = matrix(rep(0,length(subX)*2), nrow = nobs*2)
  tempX[seq(1,nobs*2,2),] = subX
  lgtdata = c(lgtdata, list(list(X = tempX,y = suby)))
  # store the data in lgtdata
  # warning c(), using list(), to prevent decompose the elements
}
mcmc <- list(R=20000,keep=1)
Z = matrix(rep(0,946))
out = rhierMnlRwMixture(Data = list(p =2,lgtdata = lgtdata, Z = Z), Prior = list(ncomp = 9), Mcmc = mcmc)

# Z here is a bug.




rhierMnlRwMixture = 
function (Data, Prior, Mcmc) 
{
  if (missing(Data)) {
    pandterm("Requires Data argument -- list of p,lgtdata, and (possibly) Z")
  }
  if (is.null(Data$p)) {
    pandterm("Requires Data element p (# choice alternatives)")
  }
  p = Data$p
  if (is.null(Data$lgtdata)) {
    pandterm("Requires Data element lgtdata (list of data for each unit)")
  }
  lgtdata = Data$lgtdata
  nlgt = length(lgtdata)
  drawdelta = TRUE
  if (is.null(Data$Z)) {
    cat("Z not specified", fill = TRUE)
    fsh()
    drawdelta = FALSE
  }
  else {
    if (!is.matrix(Data$Z)) {
      pandterm("Z must be a matrix")
    }
    else {
      if (nrow(Data$Z) != nlgt) {
        pandterm(paste("Nrow(Z) ", nrow(Z), "ne number logits ", 
                       nlgt))
      }
      else {
        Z = Data$Z
      }
    }
  }
  if (drawdelta) {
    nz = ncol(Z)
    colmeans = apply(Z, 2, mean)
    if (sum(colmeans) > 1e-05) {
      pandterm(paste("Z does not appear to be de-meaned: colmeans= ", 
                     colmeans))
    }
  }
  ypooled = NULL
  Xpooled = NULL
  if (!is.null(lgtdata[[1]]$X & is.matrix(lgtdata[[1]]$X))) {
    oldncol = ncol(lgtdata[[1]]$X)
  }
  for (i in 1:nlgt) {
    if (is.null(lgtdata[[i]]$y)) {
      pandterm(paste0("Requires element y of lgtdata[[", 
                      i, "]]"))
    }
    if (is.null(lgtdata[[i]]$X)) {
      pandterm(paste0("Requires element X of lgtdata[[", 
                      i, "]]"))
    }
    if (!is.matrix(lgtdata[[i]]$X)) {
      pandterm(paste0("lgtdata[[", i, "]]$X must be a matrix"))
    }
    if (!is.vector(lgtdata[[i]]$y, mode = "numeric") & !is.vector(lgtdata[[i]]$y, 
                                                                  mode = "logical") & !is.matrix(lgtdata[[i]]$y)) {
      pandterm(paste0("lgtdata[[", i, "]]$y must be a numeric or logical vector or matrix"))
    }
    if (is.matrix(lgtdata[[i]]$y)) {
      if (ncol(lgtdata[[i]]$y) > 1) {
        pandterm(paste0("lgtdata[[", i, "]]$y must be a vector or one-column matrix"))
      }
    }
    ypooled = c(ypooled, lgtdata[[i]]$y)
    nrowX = nrow(lgtdata[[i]]$X)
    if ((nrowX/p) != length(lgtdata[[i]]$y)) {
      pandterm(paste("nrow(X) ne p*length(yi); exception at unit", 
                     i))
    }
    newncol = ncol(lgtdata[[i]]$X)
    if (newncol != oldncol) {
      pandterm(paste("All X elements must have same # of cols; exception at unit", 
                     i))
    }
    Xpooled = rbind(Xpooled, lgtdata[[i]]$X)
    oldncol = newncol
  }
  nvar = ncol(Xpooled)
  levely = as.numeric(levels(as.factor(ypooled)))
  if (length(levely) != p) {
    pandterm(paste("y takes on ", length(levely), " values -- must be = p"))
  }
  bady = FALSE
  for (i in 1:p) {
    if (levely[i] != i) 
      bady = TRUE
  }
  cat("Table of Y values pooled over all units", fill = TRUE)
  print(table(ypooled))
  if (bady) {
    pandterm("Invalid Y")
  }
  if (missing(Prior)) {
    pandterm("Requires Prior list argument (at least ncomp)")
  }
  if (is.null(Prior$ncomp)) {
    pandterm("Requires Prior element ncomp (num of mixture components)")
  }
  else {
    ncomp = Prior$ncomp
  }
  if (is.null(Prior$SignRes)) {
    SignRes = rep(0, nvar)
  }
  else {
    SignRes = Prior$SignRes
  }
  if (length(SignRes) != nvar) {
    pandterm("The length SignRes must be equal to the dimension of X")
  }
  if (sum(!(SignRes %in% c(-1, 0, 1)) > 0)) {
    pandterm("All elements of SignRes must be equal to -1, 0, or 1")
  }
  if (is.null(Prior$mubar) & sum(abs(SignRes)) == 0) {
    mubar = matrix(rep(0, nvar), nrow = 1)
  }
  else {
    if (is.null(Prior$mubar) & sum(abs(SignRes)) > 0) {
      mubar = matrix(rep(0, nvar) + 2 * abs(SignRes), nrow = 1)
    }
    else {
      mubar = matrix(Prior$mubar, nrow = 1)
    }
  }
  if (ncol(mubar) != nvar) {
    pandterm(paste("mubar must have ncomp cols, ncol(mubar)= ", 
                   ncol(mubar)))
  }
  if (is.null(Prior$Amu) & sum(abs(SignRes)) == 0) {
    Amu = matrix(BayesmConstant.A, ncol = 1)
  }
  else {
    if (is.null(Prior$Amu) & sum(abs(SignRes)) > 0) {
      Amu = matrix(BayesmConstant.A * 10, ncol = 1)
    }
    else {
      Amu = matrix(Prior$Amu, ncol = 1)
    }
  }
  if (ncol(Amu) != 1 | nrow(Amu) != 1) {
    pandterm("Am must be a 1 x 1 array")
  }
  if (is.null(Prior$nu) & sum(abs(SignRes)) == 0) {
    nu = nvar + BayesmConstant.nuInc
  }
  else {
    if (is.null(Prior$nu) & sum(abs(SignRes)) > 0) {
      nu = nvar + BayesmConstant.nuInc + 12
    }
    else {
      nu = Prior$nu
    }
  }
  if (nu < 1) {
    pandterm("invalid nu value")
  }
  if (is.null(Prior$V) & sum(abs(SignRes)) == 0) {
    V = nu * diag(nvar)
  }
  else {
    if (is.null(Prior$V) & sum(abs(SignRes)) > 0) {
      V = nu * diag(abs(SignRes) * 0.1 + (!abs(SignRes)) * 
                      4)
    }
    else {
      V = Prior$V
    }
  }
  if (sum(dim(V) == c(nvar, nvar)) != 2) 
    pandterm("Invalid V in prior")
  if (is.null(Prior$Ad) & drawdelta) {
    Ad = BayesmConstant.A * diag(nvar * nz)
  }
  else {
    Ad = Prior$Ad
  }
  if (drawdelta) {
    if (ncol(Ad) != nvar * nz | nrow(Ad) != nvar * nz) {
      pandterm("Ad must be nvar*nz x nvar*nz")
    }
  }
  if (is.null(Prior$deltabar) & drawdelta) {
    deltabar = rep(0, nz * nvar)
  }
  else {
    deltabar = Prior$deltabar
  }
  if (drawdelta) {
    if (length(deltabar) != nz * nvar) {
      pandterm("deltabar must be of length nvar*nz")
    }
  }
  if (is.null(Prior$a)) {
    a = rep(BayesmConstant.a, ncomp)
  }
  else {
    a = Prior$a
  }
  if (length(a) != ncomp) {
    pandterm("Requires dim(a)= ncomp (no of components)")
  }
  bada = FALSE
  for (i in 1:ncomp) {
    if (a[i] < 0) 
      bada = TRUE
  }
  if (bada) 
    pandterm("invalid values in a vector")
  if (is.null(Prior$nu) & sum(abs(SignRes)) > 0) {
    nu = nvar + 15
  }
  if (is.null(Prior$Amu) & sum(abs(SignRes)) > 0) {
    Amu = matrix(0.1)
  }
  if (is.null(Prior$V) & sum(abs(SignRes)) > 0) {
    V = nu * (diag(nvar) - diag(abs(SignRes) > 0) * 0.8)
  }
  if (missing(Mcmc)) {
    pandterm("Requires Mcmc list argument")
  }
  else {
    if (is.null(Mcmc$s)) {
      s = BayesmConstant.RRScaling/sqrt(nvar)
    }
    else {
      s = Mcmc$s
    }
    if (is.null(Mcmc$w)) {
      w = BayesmConstant.w
    }
    else {
      w = Mcmc$w
    }
    if (is.null(Mcmc$keep)) {
      keep = BayesmConstant.keep
    }
    else {
      keep = Mcmc$keep
    }
    if (is.null(Mcmc$R)) {
      pandterm("Requires R argument in Mcmc list")
    }
    else {
      R = Mcmc$R
    }
    if (is.null(Mcmc$nprint)) {
      nprint = BayesmConstant.nprint
    }
    else {
      nprint = Mcmc$nprint
    }
    if (nprint < 0) {
      pandterm("nprint must be an integer greater than or equal to 0")
    }
  }
  cat(" ", fill = TRUE)
  cat("Starting MCMC Inference for Hierarchical Logit:", fill = TRUE)
  cat("   Normal Mixture with", ncomp, "components for first stage prior", 
      fill = TRUE)
  cat(paste("  ", p, " alternatives; ", nvar, " variables in X"), 
      fill = TRUE)
  cat(paste("   for ", nlgt, " cross-sectional units"), fill = TRUE)
  cat(" ", fill = TRUE)
  cat("Prior Parms: ", fill = TRUE)
  cat("nu =", nu, fill = TRUE)
  cat("V ", fill = TRUE)
  print(V)
  cat("mubar ", fill = TRUE)
  print(mubar)
  cat("Amu ", fill = TRUE)
  print(Amu)
  cat("a ", fill = TRUE)
  print(a)
  if (drawdelta) {
    cat("deltabar", fill = TRUE)
    print(deltabar)
    cat("Ad", fill = TRUE)
    print(Ad)
  }
  if (sum(abs(SignRes)) != 0) {
    cat("Sign Restrictions Vector (0: unconstrained, 1: positive, -1: negative)", 
        fill = TRUE)
    print(matrix(SignRes, ncol = 1))
  }
  cat(" ", fill = TRUE)
  cat("MCMC Parms: ", fill = TRUE)
  cat(paste("s=", round(s, 3), " w= ", w, " R= ", R, " keep= ", 
            keep, " nprint= ", nprint), fill = TRUE)
  cat("", fill = TRUE)
  oldbetas = matrix(double(nlgt * nvar), ncol = nvar)
  llmnlFract = function(beta, y, X, betapooled, rootH, w, wgt, 
                        SignRes = rep(0, ncol(X))) {
    z = as.vector(rootH %*% (beta - betapooled))
    return((1 - w) * llmnl_con(beta, y, X, SignRes) + w * 
             wgt * (-0.5 * (z %*% z)))
  }
  mnlHess_con = function(betastar, y, X, SignRes = rep(0, ncol(X))) {
    beta = betastar
    beta[SignRes != 0] = SignRes[SignRes != 0] * exp(betastar[SignRes != 
                                                                0])
    n = length(y)
    j = nrow(X)/n
    k = ncol(X)
    Xbeta = X %*% beta
    Xbeta = matrix(Xbeta, byrow = T, ncol = j)
    Xbeta = exp(Xbeta)
    iota = c(rep(1, j))
    denom = Xbeta %*% iota
    Prob = Xbeta/as.vector(denom)
    Hess = matrix(double(k * k), ncol = k)
    for (i in 1:n) {
      p = as.vector(Prob[i, ])
      A = diag(p) - outer(p, p)
      Xt = X[(j * (i - 1) + 1):(j * i), ]
      Hess = Hess + crossprod(Xt, A) %*% Xt
    }
    lambda = c(rep(1, length(SignRes)))
    lambda[SignRes == 1] = beta[SignRes == 1]
    lambda[SignRes == -1] = -beta[SignRes == -1]
    Hess = Hess * crossprod(t(lambda))
    return(Hess)
  }
  cat("initializing Metropolis candidate densities for ", nlgt, 
      " units ...", fill = TRUE)
  fsh()
  betainit = c(rep(0, nvar))
  noRes = c(rep(0, nvar))
  out = optim(betainit, llmnl_con, method = "BFGS", control = list(fnscale = -1, 
                                                                   trace = 0, reltol = 1e-06), X = Xpooled, y = ypooled, 
              SignRes = noRes)
  betainit = out$par
  betainit[SignRes != 0] = 0
  out = optim(betainit, llmnl_con, control = list(fnscale = -1, 
                                                  trace = 0, reltol = 1e-06), X = Xpooled, y = ypooled, 
              SignRes = SignRes)
  betapooled = out$par
  if (sum(abs(betapooled[SignRes]) > 10)) {
    cat("In tuning Metropolis algorithm, constrained pooled parameter estimates contain very small values", 
        fill = TRUE)
    print(cbind(betapooled, SignRes))
    cat("check any constrained values with absolute value > 10 above - implies abs(beta) > exp(10)", 
        fill = TRUE)
  }
  H = mnlHess_con(betapooled, ypooled, Xpooled, SignRes)
  rootH = chol(H)
  for (i in 1:nlgt) {
    wgt = length(lgtdata[[i]]$y)/length(ypooled)
    out = optim(betapooled, llmnlFract, method = "BFGS", 
                control = list(fnscale = -1, trace = 0, reltol = 1e-04), 
                X = lgtdata[[i]]$X, y = lgtdata[[i]]$y, betapooled = betapooled, 
                rootH = rootH, w = w, wgt = wgt, SignRes = SignRes)
    if (out$convergence == 0) {
      hess = mnlHess_con(out$par, lgtdata[[i]]$y, lgtdata[[i]]$X, 
                         SignRes)
      lgtdata[[i]] = c(lgtdata[[i]], list(converge = 1, 
                                          betafmle = out$par, hess = hess))
    }
    else {
      lgtdata[[i]] = c(lgtdata[[i]], list(converge = 0, 
                                          betafmle = c(rep(0, nvar)), hess = diag(nvar)))
    }
    oldbetas[i, ] = lgtdata[[i]]$betafmle
    if (i%%50 == 0) 
      cat("  completed unit #", i, fill = TRUE)
    fsh()
  }
  ind = NULL
  ninc = floor(nlgt/ncomp)
  for (i in 1:(ncomp - 1)) {
    ind = c(ind, rep(i, ninc))
  }
  if (ncomp != 1) {
    ind = c(ind, rep(ncomp, nlgt - length(ind)))
  }
  else {
    ind = rep(1, nlgt)
  }
  oldprob = rep(1/ncomp, ncomp)
  if (drawdelta) {
    olddelta = rep(0, nz * nvar)
  }
  else {
    olddelta = 0
    Z = matrix(0)
    deltabar = 0
    Ad = matrix(0)
  }
  draws = rhierMnlRwMixture_rcpp_loop(lgtdata, Z, deltabar, 
                                      Ad, mubar, Amu, nu, V, s, R, keep, nprint, drawdelta, 
                                      as.matrix(olddelta), a, oldprob, oldbetas, ind, SignRes)
  if (drawdelta) {
    attributes(draws$Deltadraw)$class = c("bayesm.mat", "mcmc")
    attributes(draws$Deltadraw)$mcpar = c(1, R, keep)
  }
  attributes(draws$betadraw)$class = c("bayesm.hcoef")
  attributes(draws$nmix)$class = "bayesm.nmix"
  return(draws)
}