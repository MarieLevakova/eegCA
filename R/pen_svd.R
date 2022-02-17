# Regularized SVD structure of the matrix Pi for a given cointegration rank

pen.svd <- function(X, r, dt = 1, crit = "CV", n.lambda = 101,
                    lambda.min = 1e-12, thresh = 1e-6, n.cv = 20){

  Y <- diff(X)
  N <- dim(Y)[1]
  p <- dim(Y)[2]
  Z <- X[1:N,]

  # Standardize variables
  meanY <- apply(Y, 2, mean)
  sdY <- apply(Y, 2, sd)
  Ystd <- (Y - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanY)))# %*% diag(1/sdY)

  meanZ <- apply(Z, 2, mean)
  sdZ <- apply(Z, 2, sd)
  Zstd <- (Z - matrix(1, nrow = N, ncol = 1) %*% t(as.matrix(meanZ)))# %*% diag(1/sdZ)

  fit.coint <- johansenAggregated(Y = Ystd, Z = Zstd,
                                  r = r, dt = dt, intercept = T, normalize = F)
  Pi.coint <- matrix(fit.coint$alpha, ncol = r) %*% t(matrix(fit.coint$beta, ncol = r)) * dt

  svd0 <- svd(t(Pi.coint))
  svd0daux <- matrix(0, nrow = r, ncol = r)
  diag(svd0daux) <- svd0$d[1:r]

  svd0$d <- svd0daux
  svd0$u <- matrix(svd0$u[,1:r], nrow = p, ncol = r)
  svd0$v <- matrix(svd0$v[,1:r], nrow = p, ncol = r)

  # Define weights
  w.d <- (diag(svd0$d))^(-2)
  w.u <- (svd0$u)^(-2)
  w.v <- (svd0$v)^(-2)

  # Determine the sequence of lambdas
  lambda.max <- rep(NA,r)
  ij.max <- matrix(NA, nrow = r, ncol = 2)
  sgn.max <- rep(NA, r)
  for(i in 1:r){
    # Construct the layer
    C.k <- svd0$d[i,i]*matrix(svd0$u[,i], ncol = 1)%*%matrix(svd0$v[,i], nrow = 1)
    Y.k <- Ystd - Zstd%*%(t(Pi.coint) - C.k)
    Wij <- w.d[i] * matrix(w.u[,i], ncol = 1) %*% matrix(w.v[,i], nrow = 1)
    lambda.max[i] <- max(abs((t(Y.k) %*% Zstd) * (1/Wij)))
    ij.max[i,] <- which(abs((t(Y.k) %*% Zstd) * (1/Wij)) == lambda.max[i], arr.ind = T)
    sgn.max[i] <- sign(((t(Y.k) %*% Zstd) * (1/Wij))[which.max(abs((t(Y.k) %*% Zstd) * (1/Wij)))])
  }
  lambda.seq <- 10^seq(log10(max(lambda.max)), log10(lambda.min), length = n.lambda);

  if(crit=="fixed"){
    lambda.opt[i.r] <- lambda
  } else {

    if(crit!="CV"){
      aic <- bic <- hq <- rep(NA, n.lambda)

      svd0.u <- svd0$u
      svd0.d <- svd0$d
      svd0.v <- svd0$v

      for(ii in 1:n.lambda){
        lambda <- lambda.seq[ii]
        Pi.svd <- penSvdLoop(Ystd, Zstd, svd0.u, svd0.d, svd0.v,
                             r, w.u, w.v, w.d, svd0, lambda, lambda.max, ij.max, sgn.max, thresh)
        Pi.hat <- Pi.svd$u %*% Pi.svd$d %*% t(Pi.svd$v)
        k <- sum(Pi.hat!=0)
        res <- Ystd - Zstd %*% t(Pi.hat)
        Omega.select <- (t(res) %*% res)/N
        Omega.inv <- solve(Omega.select)
        aic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k +
          sum(diag(res %*% Omega.inv %*% t(res)))
        bic[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + k*log(N) +
          sum(diag(res %*% Omega.inv %*% t(res)))
        hq[ii] <- N*p*log(2*pi) + N*log(det(Omega.select)) + 2*k*log(log(N)) +
          sum(diag(res %*% Omega.inv %*% t(res)))

        svd0.u <- Pi.svd$u
        svd0.d <- Pi.svd$d
        svd0.v <- Pi.svd$v
      }

      if (crit=="AIC") lambda.opt <- lambda.seq[which.min(aic)]
      if (crit=="BIC") lambda.opt <- lambda.seq[which.min(bic)]
      if (crit=="HQ") lambda.opt <- lambda.seq[which.min(hq)]
    }

    if(crit=="CV"){
      cv <- rep(0, n.lambda)

      # Run crossvalidation n.cv times
      for(cv_run in 1:n.cv){

        # Divide data into 5 folds
        folds <- sample(1:5, N, replace = T)
        for(i in 1:5){
          Ystd.cv <- Ystd[folds!=i,]
          Zstd.cv <- Zstd[folds!=i,]

          svd0.u <- svd0$u
          svd0.d <- svd0$d
          svd0.v <- svd0$v

          for(ii in 1:n.lambda){
            lambda <- lambda.seq[i]
            Pi.svd <- penSvdLoop(Ystd, Zstd, svd0.u, svd0.d, svd0.v,
                                 r, w.u, w.v, w.d, svd0, lambda, thresh)
            Pi.hat <- Pi.svd$u %*% diag(Pi.svd$d) %*% t(Pi.svd$v)
            res <- Ystd[folds==i,] - Zstd[folds==i,]*t(Pi.hat)
            cv[i] <- cv[i] + sum(diag(res %*% t(res)))/(N*p*n.cv)

            svd0.u <- Pi.svd$u
            svd0.d <- Pi.svd$d
            svd0.v <- Pi.svd$v
          }
        }
      }

      lambda.opt <- lambda.seq[which.min(cv)]
    }
  }

  # Fit the final model
  # (Go through the whole sequence of lambdas to achieve warm starts)
  svd0.u <- svd0$u
  svd0.d <- svd0$d
  svd0.v <- svd0$v
  for(ii in 1:n.lambda){
    lambda <- lambda.seq[i]
    Pi.svd <- penSvdLoop(Ystd, Zstd, svd0.u, svd0.d, svd0.v,
                         r, w.u, w.v, w.d, svd0, lambda, thresh)
    Pi.hat <- Pi.svd$u %*% diag(Pi.svd$d) %*% t(Pi.svd$v)
    svd0.u <- Pi.svd$u
    svd0.d <- Pi.svd$d
    svd0.v <- Pi.svd$v
    if(lambda == lambda.opt){
      Pi.final = Pi.hat
    }
  }

  mu.hat <- meanY - Pi.final %*% as.matrix(meanZ, ncol =1)
  res <- Y - matrix(1, nrow = N, ncol = 1) %*% t(mu.hat) - Z %*% t(Pi.final)
  Omega.hat <- (t(res) %*% res)/N

  return(list(PI = Pi.final/dt, MU = mu.hat/dt, OMEGA = Omega.hat/dt, lambda = lambda.opt))
}

penSvdLoop <- function(Ystd, Zstd, svd.u, svd.d, svd.v, r, w.u, w.v, w.d,
                       svd0, lambda, lambda.max, ij.max, sgn.max, thresh){

  # Exclusive extraction algorithm

  # Loop over the r layers
  for(k in 1:r){

    # Construct the layer
    C.k <- svd0$d[k,k]*matrix(svd0$u[,k], ncol = 1)%*%matrix(svd0$v[,k], nrow = 1)
    Y.k <- Ystd - Zstd%*%(svd0$u %*% svd0$d %*% t(svd0$v) - C.k)

    if(lambda < lambda.max[k]){
      # Initialize u and v
      u.new <- svd.u[,k]
      v.new <- svd.v[,k]

      # Initialize with nonzero values, if lambda.max has just been crossed
      if(all(u.new == 0) & all(v.new == 0)){
        u.new[ij.max[k,1]] <- 1
        v.new[ij.max[k,2]] <- sgn.max[k]
      }

      conv.crit <- 10*thresh
      Y4u <- matrix(Ystd, ncol = 1)
      while(conv.crit>thresh){

        # Optimize wrt v
        Z4v <- kronecker(diag(1/w.v[,k]), Zstd %*% u.new)
        lambda4v <- lambda*w.d[k]*sum(w.u[,k]*abs(svd.u[,k]))
        LASSO4v <- glmnet(y = Y4u, x = Z4v, intercept = F, family = "gaussian")
        v.hat <- as.matrix(coef(LASSO4v, s = lambda4v)[-1])
        d.hat <- sqrt(sum(diag(1/w.v[,k]) %*% v.hat)^2)
        v.new <- diag(1/(d.hat*w.v[,k])) %*% v.hat

        # Optimize wrt u
        Z4u <- kronecker(matrix(v.new, ncol = 1), Zstd) %*% diag(1/w.u[,k])
        lambda4u <- lambda*w.d[k]*sum(w.v[,k]*abs(svd.v[,k]))
        LASSO4u <- glmnet(y = Y4u, x = Z4u, intercept = F, family = "gaussian")
        u.hat <- as.matrix(coef(LASSO4u, s = lambda4u)[-1])
        d.hat <- sqrt(sum(diag(1/w.u[,k]) %*% u.hat)^2)
        u.new <- diag(1/(d.hat*w.u[,k])) %*% u.hat

        C.new <- d.hat * u.new %*% t(v.new)
        conv.crit <- sqrt(sum((C.new-C.k)^2))/sqrt(sum((C.k)^2))
        C.k <- C.new
        print(conv.crit)
      }

      svd.u[,k] <- u.new
      svd.v[,k] <- v.new
      svd.d[k,k] <- d.hat
    } else {
      svd.u[,k] <- rep(0, p)
      svd.v[,k] <- rep(0, p)
    }
  }

  return(list(d = svd.d, u = svd.u, v = svd.v))
}
