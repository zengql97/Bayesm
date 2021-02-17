myruniregGibbs=
  function (Data, Prior, Mcmc) 
  {
    pandterm = function(message) {
      stop(message, call. = FALSE)
    }
    if (missing(Data)) {
      pandterm("Requires Data argument -- list of y and X")
    }
    if (is.null(Data$X)) {
      pandterm("Requires Data element X")
    }
    X = Data$X
    if (is.null(Data$y)) {
      pandterm("Requires Data element y")
    }
    y = Data$y		
    nvar = ncol(X)   	# the number of independent variables
    nobs = length(y) 	# the number of observations
    if (nobs != nrow(X)) {
      pandterm("length(y) ne nrow(X)")
    }
    if (missing(Prior)) {
      betabar = c(rep(0, nvar))
      A = diag(rep(0.01, nvar))
      nu = 3
      ssq = var(y)
    }
    else {
      if (is.null(Prior$betabar)) {
        betabar = c(rep(0, nvar))
      }
      else {
        betabar = Prior$betabar
      }
      if (is.null(Prior$A)) {
        A = diag(rep(0.01, nvar))
      }
      else {
        A = Prior$A
      }
      if (is.null(Prior$nu)) {
        nu = 3
      }
      else {
        nu = Prior$nu
      }
      if (is.null(Prior$ssq)) {
        ssq = var(y)
      }
      else {
        ssq = Prior$ssq
      }
    }
    if (ncol(A) != nrow(A) || ncol(A) != nvar || nrow(A) != nvar) {
      pandterm(paste("bad dimensions for A", dim(A)))
    }
    if (length(betabar) != nvar) {
      pandterm(paste("betabar wrong length, length= ", length(betabar)))
    }
    if (missing(Mcmc)) {
      pandterm("requires Mcmc argument")
    }
    else {
      if (is.null(Mcmc$R)) {
        pandterm("requires Mcmc element R")
      }
      else {
        R = Mcmc$R
      }
      if (is.null(Mcmc$keep)) {
        keep = 1
      }
      else {
        keep = Mcmc$keep
      }
      if (is.null(Mcmc$sigmasq)) {
        sigmasq = var(y)
      }
      else {
        sigmasq = Mcmc$sigmasq
      }
    }
    cat(" ", fill = TRUE)
    cat("Starting Gibbs Sampler for Univariate Regression Model", 
        fill = TRUE)
    cat("  with ", nobs, " observations", fill = TRUE)
    cat(" ", fill = TRUE)
    cat("Prior Parms: ", fill = TRUE)
    cat("betabar", fill = TRUE)
    print(betabar)
    cat("A", fill = TRUE)
    print(A)
    cat("nu = ", nu, " ssq= ", ssq, fill = TRUE)
    cat(" ", fill = TRUE)
    cat("MCMC parms: ", fill = TRUE)
    cat("R= ", R, " keep= ", keep, fill = TRUE)
    cat(" ", fill = TRUE)
    sigmasqdraw = double(floor(Mcmc$R/keep)) 
    betadraw = matrix(double(floor(Mcmc$R * nvar/keep)), ncol = nvar)
    XpX = crossprod(X) 	   # t(X)*X
    Xpy = crossprod(X, y)  # t(X)*y 
    #
    # Xpx and Xpy are not varied in each iteration. 
    # thus, they can be computed outside MCMC iteration. 
    #
    sigmasq = as.vector(sigmasq)   # re-define sigmasq as a vector
    itime = proc.time()[3]
    cat("MCMC Iteration (est time to end - min) ", fill = TRUE)
    
    #
    # Different from "runireg" in simulation strategy. 
    # Generate beta, then generate sigmasq 
    #
    for (rep in 1:Mcmc$R) {
      IR = backsolve(chol(XpX/sigmasq + A), diag(nvar)) 
      # compute the inverse of the root of (XpX/sigmasq+A) 
      btilde = crossprod(t(IR)) %*% (Xpy/sigmasq + A %*% betabar)
      # Different from "runireg, the posterior mean (btilde) 
      # and the posterior variance (crossprod(t(IR))) are conditional on sigmasq    
      beta = btilde + IR %*% rnorm(nvar) 
      res = y - X %*% beta
      s = t(res) %*% res
      sigmasq = (nu * ssq + s)/rchisq(1, nu + nobs)
      # the scale parameter of the inverted chi-square draw only contain the 
      # the residual of the data. 
      sigmasq = as.vector(sigmasq)
      if (rep%%100 == 0) {
        ctime = proc.time()[3]
        timetoend = ((ctime - itime)/rep) * (R - rep)
        cat(" ", rep, " (", round(timetoend/60, 1), ")", 
            fill = TRUE)
        
      }
      if (rep%%keep == 0) {
        mkeep = rep/keep
        betadraw[mkeep, ] = bet)
a
        sigmasqdraw[mkeep] = sigmasq
      }
    }
    ctime = proc.time()[3]
    cat("  Total Time Elapsed: ", round((ctime - itime)/60, 2), 
        "\n")
    return(list(betadraw = betadraw, sigmasqdraw = sigmasqdraw))
  }

Xs =matrix(rtruncnorm(30000,0,Inf,15,20), nrow = 10000)
Xs = cbind(Xs)
Beta = c(1.23,2.34,5.14)
Ys = Xs %*% Beta + rnorm(10000,0,100)+5
lm(y~x, data = list(x= Xs, y=Ys))
myruniregGibbs(list(y=Ys,X = Xs), Mcmc = list(R=100000, keep = 10000))
