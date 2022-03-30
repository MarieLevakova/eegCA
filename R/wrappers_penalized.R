# Wrapper functions for penalized estimators of cointegration matrix.
# Output here is more user friendly and useable.

# Least-squares estimate of Pi matrix
PI.OLS <- function(Y, Z, dt = 1, intercept = T){
  N <- dim(Y)[1]
  P <- dim(Y)[2]

  if(intercept){
    m0 <- lm(Y ~ Z)
    est.coef <- m0$coefficients
    est.Omega <- (t(m0$residuals) %*% m0$residuals)/N
    out <- list(PI = t(est.coef[-1,])/dt, MU = est.coef[1,]/dt, OMEGA = est.Omega/dt)
  } else {
    m0 <- lm(Y ~ Z-1)
    est.coef <- m0$coefficients
    est.Omega <- (t(m0$residuals) %*% m0$residuals)/N
    out <- list(PI = t(est.coef)/dt, MU = est.coef[1,]/dt, OMEGA = est.Omega/dt)
  }
}

# Unrestricted penalized estimation of Pi
pen.PI.unrestricted <- function(yt, crit = c("fixed", "CV", "AIC", "BIC", "HQ"),
                                lambda, n.lambda = 100, n.cv = 10, thresh = 1e-12,
                                maxit = 1e6, dt = 1){
  Y <- diff(yt)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- yt[1:N,]

  # Standardize variables
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY))) %*% diag(1/sdY)

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ))) %*% diag(1/sdZ)

  PI.Sparse <- matrix(NA, nrow = p, ncol = p)

  # Determine the lambda sequence
  lambda.max <- 0
  lambda.min <- NA
  for(i.r in 1:p){
    penalty.factor <- rep(1, p)
    penalty.factor[i.r] <- 0
    determine_lambdasequence <- glmnet(y = Ystd[,i.r], x = Zstd,
                                       intercept = F,
                                       family = "gaussian",
                                       penalty.factor = penalty.factor)
    lambda.max <- max(c(lambda.max, determine_lambdasequence$lambda))
    lambda.min <- min(c(lambda.min, determine_lambdasequence$lambda), na.rm = T)
  }

  lambda.seq <- exp(seq(log(lambda.max), log(lambda.min), length = n.lambda))

  if(crit=="fixed"){
    lambda.opt[i.r] <- lambda
  } else {

    pi.select <- array(NA, dim = c(p, p, n.lambda))
    if(crit=="AIC"){
      aic <- rep(NA, n.lambda)
      for(i.r in 1:p){
        penalty.factor <- rep(1, p)
        penalty.factor[i.r] <- 0
        determine_lambda <- glmnet(y = Ystd[,i.r], x = Zstd,
                                   intercept = F,
                                   family = "gaussian",
                                   lambda = lambda.seq,
                                   penalty.factor = penalty.factor)
        pi.select[,i.r,] <- matrix(determine_lambda$beta, nrow = p)
      }
      for(ii in 1:n.lambda){
        N <- length(Ystd[,i.r])
        k <- sum(pi.select[,,ii]!=0)
        res <- Ystd - Zstd %*% pi.select[,,ii]
        Omega.select <- (t(res) %*% res)/N
        Omega.inv <- solve(Omega.select)
        aic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k +
          sum(diag(res %*% Omega.inv %*% t(res)))
      }
      lambda.opt <- lambda.seq[which.min(aic)]
    }

    if(crit=="BIC"){
      bic <- rep(NA, n.lambda)
      for(i.r in 1:p){
        penalty.factor <- rep(1, p)
        penalty.factor[i.r] <- 0
        determine_lambda <- glmnet(y = Ystd[,i.r], x = Zstd,
                                   intercept = F,
                                   family = "gaussian",
                                   lambda = lambda.seq,
                                   penalty.factor = penalty.factor)
        pi.select[,i.r,] <- matrix(determine_lambda$beta, nrow = p)
      }
      for(ii in 1:n.lambda){
        N <- length(Ystd[,i.r])
        k <- sum(pi.select[,,ii]!=0)
        res <- Ystd - Zstd %*% pi.select[,,ii]
        Omega.select <- (t(res) %*% res)/N
        Omega.inv <- solve(Omega.select)
        bic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + k*log(N*p) +
          sum(diag(res %*% Omega.inv %*% t(res)))
      }
      lambda.opt <- lambda.seq[which.min(bic)]
    }

    if(crit=="HQ"){
      hq <- rep(NA, n.lambda)
      for(i.r in 1:p){
        penalty.factor <- rep(1, p)
        penalty.factor[i.r] <- 0
        determine_lambda <- glmnet(y = Ystd[,i.r], x = Zstd,
                                   intercept = F,
                                   family = "gaussian",
                                   lambda = lambda.seq,
                                   penalty.factor = penalty.factor)
        pi.select[,i.r,] <- matrix(determine_lambda$beta, nrow = p)
      }
      for(ii in 1:n.lambda){
        N <- length(Ystd[,i.r])
        k <- sum(pi.select[,,ii]!=0)
        res <- Ystd - Zstd %*% pi.select[,,ii]
        Omega.select <- (t(res) %*% res)/N
        Omega.inv <- solve(Omega.select)
        hq[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k*log(log(N*p)) +
          sum(diag(res %*% Omega.inv %*% t(res)))
      }
      lambda.opt <- lambda.seq[which.min(hq)]
    }

    if(crit=="CV"){
      cv <- rep(0, n.lambda)
      for(i.cv in 1:n.cv){
        for(i.r in 1:p){
          penalty.factor <- rep(1, p)
          penalty.factor[i.r] <- 0
          determine_lambda <- cv.glmnet(y = Ystd[,i.r], x = Zstd,
                                        intercept = F,
                                        family = "gaussian",
                                        lambda = lambda.seq,
                                        penalty.factor = penalty.factor)
          cv <- cv + 1/(n.cv*p)*determine_lambda$cvm
        }
      }
      lambda.opt <- lambda.seq[which.min(cv)]
    }
  }

  # Fit the final model

  for(i.r in 1:p){
    penalty.factor <- rep(1, p)
    penalty.factor[i.r] <- 0
    LASSOfinal <- glmnet(y = Ystd[,i.r], x = Zstd,
                         standardize = F,
                         intercept = F, lambda = lambda.seq,
                         family = "gaussian",
                         penalty.factor = penalty.factor)
    PI.Sparse[i.r,] <- matrix(coef(LASSOfinal, s = lambda.opt), nrow = 1)[-1] %*% diag(sdY[i.r]/sdZ)
  }

  mu.hat <- meanY - PI.Sparse %*% as.matrix(meanZ, ncol =1)
  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.hat) - Z %*% t(PI.Sparse)
  Omega.hat <- (t(res) %*% res)/N

  return(list(PI = PI.Sparse/dt, MU = mu.hat/dt, OMEGA = Omega.hat/dt, lambda = lambda.opt))
}

pen.PI.restricted <- function(X, n.lambda, lambda.min = 1e-12, r, maxiter = 100,
                              dt = 1, crit = "CV", w.auto = TRUE, n.cv = 100, q = 2){

  p <- dim(X)[2]

  # Translating crit from a string to an integer value
  crit.vector <- c("CV", "AIC", "BIC", "HQ")
  crit.int <- which(crit.vector ==  crit) - 1

  # Calling cpp function
  out <- penPiCpp(X, n.lambda, lambda.min, r, maxiter, crit.int, dt, w.auto, n.cv, q)

  # Extracting results from the cpp output
  PI <- out[1:p,1:p]
  MU <- out[1:p,p+1]
  OMEGA <- out[1:p,(p+2):(2*p+1)]
  lambda <- out[1,2*p+2]
  PI.iter <- out[(p+1):(p+maxiter-1),1:(p^2)]
  objective.iter <- out[(p+1):(p+maxiter-1),p^2+1]
  w <- out[(p+1):(p+maxiter),p^2+2]
  lambda.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),1]
  crit.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),2]
  Pi.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),3:(p^2+2)]

  return(list(PI = PI, MU = MU, OMEGA = OMEGA, lambda = lambda,
              PI.iter = PI.iter, objective.iter = objective.iter, w = w,
              lambda.seq = rev(lambda.seq), crit = crit, crit.seq = rev(crit.seq),
              Pi.seq = Pi.seq))
}

pen.PI.nuclear <- function(X, n.lambda, lambda.min = 1e-12, miniter = 10, maxiter = 100,
                              dt = 1, crit = "CV", n.cv = 100, thresh = 1e-6){

  p <- dim(X)[2]

  # Translating crit from a string to an integer value
  crit.vector <- c("CV", "AIC", "BIC", "HQ")
  crit.int <- which(crit.vector ==  crit) - 1

  # Calling cpp function
  out <- penNuclearCpp(X, n.lambda, lambda.min, miniter, maxiter, crit.int, dt, n.cv, thresh)

  # Extracting results from the cpp output
  PI <- out[1:p,1:p]
  MU <- out[1:p,p+1]
  OMEGA <- out[1:p,(p+2):(2*p+1)]
  lambda <- out[1,2*p+2]
  PI.iter <- out[(p+1):(p+maxiter-1),1:(p^2)]
  objective.iter <- out[(p+1):(p+maxiter-1),p^2+1]
  lambda.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),1]
  crit.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),2]
  Pi.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),3:(p^2+2)]
  iter.seq <- out[(p+maxiter+1):(p+maxiter+n.lambda),p^2+3]

  return(list(PI = PI, MU = MU, OMEGA = OMEGA, lambda = lambda,
              PI.iter = PI.iter, objective.iter = objective.iter,
              lambda.seq = rev(lambda.seq), crit = crit, crit.seq = rev(crit.seq),
              Pi.seq = Pi.seq, iter.seq = iter.seq))
}

# Penalty on beta
pen.beta <- function(X, r, max.iter = 10, conv = 10^-2, nlambda = 100, n.cv = 20,
                     glmnetthresh = 1e-04, dt = 1, equal.penalty = F){

  ## INPUT
  # X: vector time Series
  # r: cointegration rank
  # max.iter: maximum number of iterations
  # conv: convergence parameter
  # nlambda: number of penalty values to check
  # n.cv: number of repetitions of the crossvalidation procedure
  # glmnetthresh: tolerance parameter glmnet function
  # dt: timestep
  # equal.penalty: logical value, whether the penalty should be constant over the r cointegration vectors, or not

  ## OUTPUT
  # BETA: estimate of cointegrating vectors
  # ALPHA: estimate of adjustment coefficients
  # OMEGA: estimate of covariance matrix
  # PI: estimate of Pi
  # MU: estimate of the intercept
  # it: number of iterations
  # value.obj: value of the objective function at each iteration
  # Pi.it: estimate of Pi at each iteration

  ## START CODE
  Y <- diff(X)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- X[1:N,]

  # Centering variables, to remove the effect of intercept
  Ystd <- Y
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY)))

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ)))

  # Starting value
  # fit0 <- johansenAggregated(Y = Ystd, Z = Zstd, r = r, dt = dt, intercept = F, normalize = F)
  # beta.init <- fit0$beta
  mm <- diag(as.vector(1/sdZ^2)) %*% cov(Zstd, Ystd) %*% diag(as.vector(1/sdY^2)) %*% cov(Ystd, Zstd)
  beta.init <- matrix(eigen(mm)$vectors[,1:r], ncol = r)
  alpha.init <- t(coef(lm(Ystd ~ I(Zstd %*% beta.init)-1)))
  Pi.init <- matrix(alpha.init, ncol = r) %*% t(matrix(beta.init, ncol = r))
  mu.init <- meanY - Pi.init %*% as.matrix(meanZ, ncol = 1)
  RESID <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.init) - Z %*% t(Pi.init)
  Omega.init <- (t(RESID) %*% RESID)/N

  # Convergence parameters: initialization
  it <- 1
  diff.obj <- 10*conv
  value.obj <- matrix(NA, ncol = 1, nrow = max.iter+1)
  Pi.it <- matrix(NA, ncol = p*r, nrow = max.iter + 1)

  RESID <- Y - Z %*% t(Pi.init) - matrix(1, nrow = N, ncol = 1) %*% matrix(mu.init, nrow = 1)
  value.obj[1,] <- (1/N)*sum(diag(RESID %*% solve(Omega.init) %*% t(RESID))) + log(det(Omega.init))
  Pi.it[1,] <- matrix(alpha.init, nrow = 1)

  while((it < max.iter) &  (diff.obj > conv)){
    # Obtain Alpha
    FIT2 <- NTS.ALPHA.Procrusted(Y = Ystd, Z = Zstd,
                                 r = r, Omega = Omega.init, beta = beta.init)
    # Obtain Beta and Omega
    FIT3 <- NTS.BETA(Y = Ystd, Z = Zstd, r = r, Omega = Omega.init,
                     P = FIT2$P, alpha = FIT2$ALPHA, alphastar = FIT2$ALPHAstar,
                     nlambda = nlambda, n.cv = n.cv, glmnetthresh = glmnetthresh,
                     equal.penalty = equal.penalty)

    # Check convergence
    beta.new <- FIT3$BETA
    beta.init <- matrix(beta.new, nrow = p, ncol = r)
    alpha.init <- FIT2$ALPHA
    Pi.init <- alpha.init%*%t(beta.init)
    mu.init <- meanY - Pi.init %*% as.matrix(meanZ, ncol = 1)
    RESID <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.init) - Z %*% t(Pi.init)
    Omega.init <- (t(RESID) %*% RESID)/N

    value.obj[1+it,] <- (1/N)*sum(diag((RESID) %*% solve(Omega.init) %*% t(RESID))) + log(det(Omega.init))
    diff.obj <- abs(value.obj[1+it,] - value.obj[it,])/abs(value.obj[it,])
    Pi.it[1+it,] <- matrix(alpha.init, nrow = 1)
    it <- it + 1

  }

  out <- list(BETA = beta.init, ALPHA = alpha.init/dt, OMEGA = Omega.init/dt,
              PI = Pi.init/dt, MU = mu.init/dt, it = it, value.obj = value.obj,
              Pi.it = Pi.it)

  return(out)
}

# Penalty on alpha
pen.alpha <- function(X, r, dt = 1, equal.penalty = F, n.lambda = 100, n.cv = 20){

  ## INPUT
  # X: vector time Series
  # r: cointegration rank
  # max.iter: maximum number of iterations
  # conv: convergence parameter
  # nlambda: number of penalty values to check
  # glmnetthresh: tolerance parameter glmnet function
  # dt: timestep
  # equal.penalty: logical value, whether the penalty should be constant over the r cointegration vectors, or not

  ## OUTPUT
  # BETA: estimate of cointegrating vectors
  # ALPHA: estimate of adjustment coefficients
  # OMEGA: estimate of covariance matrix
  # PI: estimate of Pi
  # MU: estimate of the intercept
  # it: number of iterations
  # value.obj: value of the objective function at each iteration
  # Pi.it: estimate of Pi at each iteration

  ## START CODE
  Y <- diff(X)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- X[1:N,]

  # Centering variables, to remove the effect of intercept
  Ystd <- Y
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY)))

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ)))

  # Fit the usual cointegration model
  fit0 <- johansenAggregated(Y = Ystd, Z = Zstd, r = r, dt = dt, intercept = F, normalize = F)

  # Extract beta and create new set of predictors
  BETA <- matrix(fit0$beta, nrow = p, ncol = r)
  Z.r <- matrix(Zstd %*% BETA, nrow = N, ncol = r)

  # Estimate alpha with LASSO penalty

  ALPHA.Sparse <- matrix(NA, nrow = p, ncol = r)

  if(r == 1){
    ALPHA.Sparse  <- matrix(coef(lm(Ystd ~ Z.r - 1)), ncol = 1)
  } else {

    if(equal.penalty){
      # Determine the lambda sequence
      lambda.max <- 0
      lambda.min <- NA
      for(i.r in 1:p){
        determine_lambdasequence <- glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                                           family = "gaussian")
        lambda.max <- max(c(lambda.max, determine_lambdasequence$lambda))
        lambda.min <- min(c(lambda.min, determine_lambdasequence$lambda), na.rm = T)
      }

      lambda.seq <- exp(seq(log(lambda.max), log(lambda.min), length = n.lambda))
      cv <- rep(0, n.lambda)
      for(i.cv in 1:n.cv){
        for(i.r in 1:p){
          determine_lambda <- cv.glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                                        family = "gaussian", lambda = lambda.seq)
          cv <- cv + 1/(n.cv*r)*determine_lambda$cvm
        }
      }
      lambda.opt <- lambda.seq[which.min(cv)]

      # Fit the final model

      for(i.r in 1:p){
        LASSOfinal <- glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                             lambda = lambda.seq, family = "gaussian")
        ALPHA.Sparse[i.r,] <- matrix(coef(LASSOfinal, s = lambda.opt), nrow = 1)[-1]
      }
    } else {
      # Determine the lambda sequence
      for(i.r in 1:p){

        cv <- rep(0, n.lambda)
        for(i.cv in 1:n.cv){
          determine_lambda <- cv.glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                                        family = "gaussian", nlambda = n.lambda)
          cv <- cv + 1/(n.cv*r)*determine_lambda$cvm
        }
        lambda.opt <- determine_lambda$lambda[which.min(cv)]

        # Fit the final model
        LASSOfinal <- glmnet(y = Ystd[,i.r], x = Z.r, intercept = F,
                             nlambda = n.lambda, family = "gaussian")
        ALPHA.Sparse[i.r,] <- matrix(coef(LASSOfinal, s = lambda.opt), nrow = 1)[-1]
      }
    }
  }

  PI <- ALPHA.Sparse %*% t(BETA)
  MU <- meanY - PI %*% as.matrix(meanZ, ncol =1)
  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(MU) - Z %*% t(PI)
  OMEGA <- (t(res) %*% res)/N

  return(list(ALPHA = ALPHA.Sparse/dt, BETA = BETA, PI = PI/dt, OMEGA = OMEGA/dt, MU = MU/dt))
}

# Row sparse penalization of alpha
pen.alpha.rowSparse <- function(X, crit = c("fixed", "CV", "AIC", "BIC", "HQ"),
                   lambda, n.lambda = 100, n.cv = 10, thresh = 1e-12,
                   maxit = 1e6, dt = 1, psi = 2){
  Y <- diff(X)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- X[1:N,]

  # Standardize variables
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY))) #%*% diag(1/sdY)

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ))) #%*% diag(1/sdZ)

  # Get the cointegration matrix by Johansen procedure
  fit.johansen <- johansenAggregated(Ystd, Zstd, r = p, dt = dt, intercept = F, normalize = F)$PI

  # Construct a new predictor matrix
  beta.full <- fit.johansen$beta
  Z.new <- Zstd %*% beta.full

  alpha.Sparse <- matrix(NA, nrow = p, ncol = p)

  # Determine the lambda sequence
  determine_lambdasequence <- glmnet(y = Ystd, x = Z.new,
                                     intercept = F,
                                     family = "mgaussian")

  lambda.seq <- determine_lambdasequence$lambda
  n.lambda <- length(lambda.seq)

  if(crit=="fixed"){
    lambda.opt[i.r] <- lambda
  } else {


    if(crit!="CV"){
      alpha.select <- array(0, dim = c(p, p, n.lambda))
      for(i in 1:p){
        alpha.select[,i,] <- as.matrix(determine_lambdasequence$beta[[i]])
      }

      aic <- rep(NA, n.lambda)
      bic <- rep(NA, n.lambda)
      hq <- rep(NA, n.lambda)
      for(ii in 1:n.lambda){
        N <- dim(Ystd)[1]
        k <- sum((Smatalpha.select[,,ii] %*% t(beta.full))!=0)
        res <- Ystd - Z.new %*% alpha.select[,,ii]
        Omega.select <- (t(res) %*% res)/N
        Omega.inv <- solve(Omega.select)
        aic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k +
          sum(diag(res %*% Omega.inv %*% t(res)))
        bic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + k*log(N*p) +
          sum(diag(res %*% Omega.inv %*% t(res)))
        hq[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k*log(log(N*p)) +
          sum(diag(res %*% Omega.inv %*% t(res)))
      }
      if(crit=="AIC"){
        lambda.opt <- lambda.seq[which.min(aic)]
      }
      if(crit=="BIC"){
        lambda.opt <- lambda.seq[which.min(bic)]
      }
      if(crit=="HQ"){
        lambda.opt <- lambda.seq[which.min(hq)]
      }
    }
    if(crit=="CV"){
      cv <- rep(0, n.lambda)
      for(i.cv in 1:n.cv){
        determine_lambda <- cv.glmnet(y = Ystd, x = Z.new,
                                      intercept = F,
                                      family = "mgaussian")
        cv <- cv + 1/(n.cv*p)*determine_lambda$cvm
      }
      lambda.opt <- lambda.seq[which.min(cv)]
    }
  }

  # Fit the final model
  LASSOfinal <- glmnet(y = Ystd, x = Z.new,
                       intercept = F,
                       family = "mgaussian")
  alpha.Sparse <- matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    alpha.Sparse[,i] <- as.matrix(coef(LASSOfinal, s = lambda.opt)[[i]])[-1]
  }
  PI.Sparse <- alpha.Sparse %*% t(beta.full)
  # for(i in 1:p){
  #   PI.Sparse[,i] <- PI.Sparse[,i] %*% diag(sdY[i]/sdZ)
  # }

  mu.hat <- meanY - PI.Sparse %*% as.matrix(meanZ, ncol =1)
  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.hat) - Z %*% t(PI.Sparse)
  Omega.hat <- (t(res) %*% res)/N

  return(list(PI = PI.Sparse/dt, MU = mu.hat/dt, OMEGA = Omega.hat/dt, lambda = lambda.opt))
}


# Unrestricted penalized estimation of Pi
pen.qr <- function(X, crit = c("fixed", "CV", "AIC", "BIC", "HQ"),
                                lambda, n.lambda = 100, n.cv = 10, thresh = 1e-12,
                                maxit = 1e6, dt = 1, psi = 2){
  Y <- diff(X)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- X[1:N,]

  # Standardize variables
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY))) #%*% diag(1/sdY)

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ))) #%*% diag(1/sdZ)

  # Preliminary OLS estimate
  Pi0 <- PI.OLS(Ystd, Zstd, dt = dt, intercept = F)$PI
  QR.Pi0 <- qr(t(Pi0))
  Smat <- qr.Q(QR.Pi0)
  Rmat <- qr.R(QR.Pi0)

  # Construct a new predictor matrix
  ZS <- Zstd %*%Smat

  R.Sparse <- matrix(NA, nrow = p, ncol = p)

  # Determine the lambda sequence
  penalty.factor <- rep(NA, p)
  for(i in 1:p){
    penalty.factor[i] <- sqrt(sum((Rmat[i,i:p])^2))
  }
  determine_lambdasequence <- glmnet(y = Ystd, x = ZS,
                                     intercept = F,
                                     penalty.factor = penalty.factor^psi,
                                     family = "mgaussian")

  lambda.seq <- determine_lambdasequence$lambda
  n.lambda <- length(lambda.seq)

  if(crit=="fixed"){
    lambda.opt[i.r] <- lambda
  } else {


    if(crit!="CV"){
      R.select <- array(0, dim = c(p, p, n.lambda))
      for(i in 1:p){
        R.select[,i,] <- as.matrix(determine_lambdasequence$beta[[i]])
      }

      aic <- rep(NA, n.lambda)
      bic <- rep(NA, n.lambda)
      hq <- rep(NA, n.lambda)
      for(ii in 1:n.lambda){
        N <- dim(Ystd)[1]
        k <- sum((Smat %*% R.select[,,ii])!=0)
        res <- Ystd - ZS %*% R.select[,,ii]
        Omega.select <- (t(res) %*% res)/N
        Omega.inv <- solve(Omega.select)
        aic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k +
          sum(diag(res %*% Omega.inv %*% t(res)))
        bic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + k*log(N*p) +
          sum(diag(res %*% Omega.inv %*% t(res)))
        hq[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k*log(log(N*p)) +
          sum(diag(res %*% Omega.inv %*% t(res)))
      }
      if(crit=="AIC"){
        lambda.opt <- lambda.seq[which.min(aic)]
      }
      if(crit=="BIC"){
        lambda.opt <- lambda.seq[which.min(bic)]
      }
      if(crit=="HQ"){
        lambda.opt <- lambda.seq[which.min(hq)]
      }
    }
    if(crit=="CV"){
      cv <- rep(0, n.lambda)
      for(i.cv in 1:n.cv){
        determine_lambda <- cv.glmnet(y = Ystd, x = ZS,
                                      intercept = F,
                                      penalty.factor = penalty.factor^psi,
                                      family = "mgaussian")
        cv <- cv + 1/(n.cv*p)*determine_lambda$cvm
      }
      lambda.opt <- lambda.seq[which.min(cv)]
    }
  }

  # Fit the final model
  LASSOfinal <- glmnet(y = Ystd, x = ZS,
                       intercept = F,
                       penalty.factor = penalty.factor^psi,
                       family = "mgaussian")
  R.Sparse <- matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    R.Sparse[,i] <- as.matrix(coef(LASSOfinal, s = lambda.opt)[[i]])[-1]
  }
  PI.Sparse <- t(Smat %*% R.Sparse)
  # for(i in 1:p){
  #   PI.Sparse[,i] <- PI.Sparse[,i] %*% diag(sdY[i]/sdZ)
  # }

  mu.hat <- meanY - PI.Sparse %*% as.matrix(meanZ, ncol =1)
  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.hat) - Z %*% t(PI.Sparse)
  Omega.hat <- (t(res) %*% res)/N

  return(list(PI = PI.Sparse/dt, MU = mu.hat/dt, OMEGA = Omega.hat/dt, lambda = lambda.opt))
}

pen.PI.rank <- function(X, n.lambda, lambda.min = 1e-12, dt = 1, crit = "CV", n.cv = 100){

  p <- dim(X)[2]

  # Translating crit from a string to an integer value
  crit.vector <- c("CV", "AIC", "BIC", "HQ")
  crit.int <- which(crit.vector ==  crit) - 1

  # Calling cpp function
  out <- penRankCpp(X, n.lambda, lambda.min, crit.int, dt, n.cv)

  # Extracting results from the cpp output
  PI <- out[1:p,1:p]
  MU <- out[1:p,p+1]
  OMEGA <- out[1:p,(p+2):(2*p+1)]
  lambda <- out[1,2*p+2]
  lambda.seq <- out[(p+1):(p+n.lambda), 1]
  crit.seq <- out[(p+1):(p+n.lambda), 2]
  PI.seq <- out[(p+n.lambda):(p+1), 3:(p^2+2)] # reverse the order

  return(list(PI = PI, MU = MU, OMEGA = OMEGA, lambda = lambda,
              lambda.seq = rev(lambda.seq), crit = crit, crit.seq = rev(crit.seq),
              PI.seq = PI.seq))
}

pen.PI.nuclearAdapt <- function(X, n.lambda, lambda.min = 1e-12, dt = 1, crit = "CV",
                                n.cv = 100, w.gamma = 0){

  p <- dim(X)[2]

  # Translating crit from a string to an integer value
  crit.vector <- c("CV", "AIC", "BIC", "HQ")
  crit.int <- which(crit.vector ==  crit) - 1

  # Calling cpp function
  out <- penAdaptNuclearCpp(X, n.lambda, lambda.min, crit.int, dt, n.cv, w.gamma)

  # Extracting results from the cpp output
  PI <- out[1:p,1:p]
  MU <- out[1:p,p+1]
  OMEGA <- out[1:p,(p+2):(2*p+1)]
  lambda <- out[1,2*p+2]
  lambda.seq <- out[(p+1):(p+n.lambda), 1]
  crit.seq <- out[(p+1):(p+n.lambda), 2]

  return(list(PI = PI, MU = MU, OMEGA = OMEGA, lambda = lambda,
              lambda.seq = rev(lambda.seq), crit = crit, crit.seq = rev(crit.seq)))
}
