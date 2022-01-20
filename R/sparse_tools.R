# Auxiliary functions for solving the problem of penalized beta
NTS.ALPHA.Procrusted<-function(Y, Z, r, Omega, beta){
  ### FUNCTION TO ESTIMATE ALPHA ###

  ## INPUT
  # Y: Response Time Series
  # Z: Predictor time Series in Levels
  # r: cointegration rank
  # Omega: estimate of error covariance matrix
  # beta: estimate of cointegrating vector

  ## OUTPUT
  # ALPHA: estimate of adjustment coefficients
  # ALPHAstar: estimate of transformed adjustment coefficients
  # P: transformation matrix

  ## START CODE
  # Data matrices
  Xmatrix <- Z%*%beta
  decomp <- eigen(Omega)
  P <- (decomp$vectors) %*% diag(1/sqrt(decomp$values)) %*% solve(decomp$vectors)
  Pmin <- (decomp$vectors) %*% diag(sqrt(decomp$values)) %*% solve(decomp$vectors)

  # Singular Value Decomposition to compute ALPHA
  SingDecomp <- svd(t(Xmatrix) %*% Y %*% P)
  ALPHAstar <- t(SingDecomp$u %*% t(SingDecomp$v))
  if (r==1){
    ALPHAstar <- matrix(ALPHAstar, ncol = 1)
  }
  ALPHA <- Pmin %*% ALPHAstar

  out <- list(ALPHA = ALPHA, ALPHAstar = ALPHAstar, P = P)
  return(out)
}

NTS.BETA <- function(Y, Z, r, Omega, P, alpha, alphastar,
                     nlambda, n.cv, glmnetthresh = 1e-04, equal.penalty){
  ### FUNCTIONS TO ESTIMATE BETAs ###
  # First step: Determine Beta conditional on Gamma, alpha and Omega using Lasso Regression
  # Second step: Determine Omega conditional on Gamma, alpha and beta using GLasso

  ## INPUT
  # Y: Response Time Series
  # Z: Time Series in Levels
  # r: cointegration rank
  # Omega: estimate of error covariance matrix
  # P: transformation matrix P derived from Omega
  # alpha: estimate of adjustment coefficients
  # alphastar: estimate of transformed adjustment coefficients
  # nlambda: number of tuning parameter values to try
  # rho.glasso: tuning parameter inverse error covariance matrix
  # cutoff: cutoff value time series cross-validation approach
  # glmnetthresh: tolerance parameter glmnet function

  ## OUTPUT
  # BETA: estimate of cointegrating vectors

  # Data matrices
  Ymatrix <- Y %*% t(P) %*% alphastar
  Xmatrix <- Z
  n <- nrow(Ymatrix)

  # Store Estimates of cointegrating vector
  BETA.Sparse <- matrix(NA, ncol = r, nrow = ncol(Y))

  # Standardized Response
  Ymatrixsd <- Ymatrix
  for (i.r in 1:r){
    Ymatrixsd[,i.r] <- Ymatrix[,i.r]/sd(Ymatrix[,i.r])
  }

  if(equal.penalty){

    # Determine the lambda sequence
    lambda.max <- 0
    lambda.min <- NA
    for(i.r in 1:r){
      determine_lambdasequence <- glmnet(y = Ymatrixsd[,i.r], x = Xmatrix,
                                         intercept = F,
                                         family = "gaussian")
      lambda.max <- max(c(lambda.max, determine_lambdasequence$lambda))
      lambda.min <- min(c(lambda.min, determine_lambdasequence$lambda), na.rm = T)
    }

    lambda_restricted <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))

    cv <- rep(0, nlambda)
    for(i in 1:n.cv){
      for (i.r in 1:r){
        # Calculate crossvalidated error
        BETA.scaled <- cv.glmnet(y = Ymatrixsd[,i.r], x = Xmatrix, lambda = lambda_restricted,
                                 standardize = F, intercept = F, family = "gaussian",
                                 thresh = glmnetthresh)
        cv <- cv + 1/(n.cv*r)*BETA.scaled$cvm
      }
    }

    # Lasso with fixed lambda
    lambda.opt <- lambda_restricted[which.min(cv)]

    for(i.r in 1:r){
      LASSOfinal <- glmnet(y = Ymatrixsd[,i.r], x = Xmatrix, standardize = F, intercept = F,
                           lambda = lambda_restricted, family = "gaussian", thresh = glmnetthresh)
      BETAscaled <- matrix(coef(LASSOfinal, s = lambda.opt), nrow = 1)[-1]
      BETA.Sparse[,i.r] <- BETAscaled*sd(Ymatrix[,i.r])
    }
  } else {

    for(i.r in 1:r){
    # Determine lambda
      cv <- rep(0, nlambda)
      for(i in 1:n.cv){
        BETA.scaled <- cv.glmnet(y = Ymatrixsd[,i.r], x = Xmatrix, nlambda = nlambda,
                               standardize = F, intercept = F, family = "gaussian",
                               thresh = glmnetthresh)
        cv <- cv + 1/(n.cv*r)*BETA.scaled$cvm
      }
      # LASSO with varying lambdas
      lambda.opt <- BETA.scaled$lambda[which.min(cv)]
      LASSOfinal <- glmnet(y = Ymatrixsd[,i.r], x = Xmatrix, standardize = F, intercept = F,
                           nlambda = nlambda, family = "gaussian", thresh = glmnetthresh)
      BETAscaled <- matrix(coef(LASSOfinal, s = lambda.opt), nrow = 1)[-1]
      BETA.Sparse[,i.r] <- BETAscaled*sd(Ymatrix[,i.r])
    }
  }

  out <- list(BETA = BETA.Sparse)
}

# huge.glasso.eegCA <- function (x, lambda = NULL, cov.output = FALSE) {
#   n = nrow(x)
#   d = ncol(x)
#   x = scale(x)
#   S = cor(x)
#   rm(x)
#
#   scr <- FALSE
#   if (is.null(lambda)) {
#     nlambda = 10
#     lambda.min.ratio = 0.1
#     lambda.max = max(max(S - diag(d)), -min(S - diag(d)))
#     lambda.min = lambda.min.ratio * lambda.max
#     lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
#   }
#   fit = .Call("_huge_hugeglasso", S, lambda, scr, verbose,
#               cov.output)
#   fit$scr = scr
#   fit$lambda = lambda
#   fit$cov.input = cov.input
#   fit$cov.output = cov.output
#   rm(S)
#   return(fit)
# }

huge.eegCA <- function (x, lambda = NULL, scr.num = NULL, sym = "or") {
  gcinfo(FALSE)
  est = list()
  est$method = "glasso"
  fit = huge.glasso(x, lambda = lambda, cov.output = T, verbose = F)
  est$path = fit$path
  est$lambda = fit$lambda
  est$icov = fit$icov
  est$df = fit$df
  est$sparsity = fit$sparsity
  est$loglik = fit$loglik
  est$cov = fit$cov
  est$cov.input = fit$cov.input
  est$cov.output = fit$cov.output
  est$scr = fit$scr
  rm(fit)
  gc()
  est$data = x
  gc()
  class(est) = "huge"
  return(est)
}

huge.select.eegCA <- function (est, ebic.gamma = 0.5, rep.num = 20){
  gcinfo(FALSE)
  criterion = "ebic"
  n = nrow(est$data)
  d = ncol(est$data)
  nlambda = length(est$lambda)
  est$ebic.score = -n * est$loglik + log(n) * est$df +
    4 * ebic.gamma * log(d) * est$df
  est$opt.index = which.min(est$ebic.score)
  est$refit = est$path[[est$opt.index]]
  est$opt.icov = est$icov[[est$opt.index]]
  est$opt.cov = est$cov[[est$opt.index]]
  est$opt.lambda = est$lambda[est$opt.index]
  est$opt.sparsity = est$sparsity[est$opt.index]
  est$criterion = criterion
  class(est) = "select"
  return(est)
}
