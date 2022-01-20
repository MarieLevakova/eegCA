#####################################################
#
#  FUNCTIONS FOR DETERMINING THE COINTEGRATION RANK
#
#####################################################

rank.johansen <- function(X, conf.level = 0.05, type = c("trace", "max")){
  # A function estimating the cointegration rank by Johansen procedure,
  # namely the trace and the maximum eigenvalue test of the rank
  # Input: X - matrix with signals in columns
  #        conf.level - confidence level for the tests
  #        type - type of the test (trace test, maximum eigenvalue test)
  # Output: list with elements
  #         $r - one- or two-element vector with estimated rank(s)
  #         $lambdas - eigenvalues of the eigenvalue problem solved in the Johansen procedure
  #         $lrts - test statistics for testing all possible ranks

  alpha.levels <- c(0.1, 0.05, 0.01)
  P <- dim(X)[2]  # no. of columns

  if (P > 10) print("Dimension exceed the allowed maximum of 10.") else {
    N <- dim(X)[1]  # no. of rows

    cvals <- array(c(6.5, 12.91, 18.9, 24.78, 30.84, 36.25,
                     42.06, 48.43, 54.01, 59, 65.07, 8.18, 14.9, 21.07, 27.14,
                     33.32, 39.43, 44.91, 51.07, 57, 62.42, 68.27, 11.65,
                     19.19, 25.75, 32.14, 38.78, 44.59, 51.3, 57.07, 63.37,
                     68.61, 74.36, 6.5, 15.66, 28.71, 45.23, 66.49, 85.18,
                     118.99, 151.38, 186.54, 226.34, 269.53, 8.18, 17.95,
                     31.52, 48.28, 70.6, 90.39, 124.25, 157.11, 192.84, 232.49,
                     277.39, 11.65, 23.52, 37.22, 55.43, 78.87, 104.2, 136.06,
                     168.92, 204.79, 246.27, 292.65), c(11, 3, 2))

    lambdas <- vecmEigen(X)
    lrts1 <- -N*rev(cumsum(log(1-rev(lambdas))))
    lrts2 <- -N*log(1-lambdas)

    r <- c(NA, NA)
    names(r) <- c("trace", "max")
    r["trace"] <- which(lrts1 <= cvals[0:P, conf.level == alpha.levels, 2])[1] - 1
    r["max"] <- which(lrts2 <= cvals[0:P, conf.level == alpha.levels, 1])[1] - 1
    return(list(r = r[type],
         lambdas = lambdas,
         lrts = data.frame(trace = lrts1, max = lrts2, row.names = paste("r =", 0:(P-1)))))
  }
}


# Johansen procedure with bootstrap method
rank.bootstrap <- function(X, B, dt, conf.level = 0.05){
  bboot <- suppressWarnings(bootstrap(X, B, dt))
  lower <- apply(bboot$boot, 2, quantile, probs = conf.level/2, na.rm = T)
  upper <- apply(bboot$boot, 2, quantile, probs = 1-conf.level/2, na.rm = T)
  which(bboot$test <= upper)[1] - 1
}

rank.bootstrapAggregated <- function(Y, Z, B, dt, conf.level = 0.05, n=nrow(Y), n.epoch=1){
  bboot <- suppressWarnings(bootstrapAggregated(Y, Z, B, dt, n = n, n.epoch = n.epoch))
  lower <- apply(bboot$boot, 2, quantile, probs = conf.level/2)
  upper <- apply(bboot$boot, 2, quantile, probs = 1-conf.level/2)
  which(bboot$test <= upper)[1] - 1
}

# Method by Bunea et al.
rank.bunea <- function(X, figure = FALSE){
  p <- ncol(X)
  N <- nrow(X)
  l <- rankMatrix(X)[1]

  Delta_phi <- diff(X, differences = 1)
  R0 <- lm(Delta_phi ~ 1)$residuals
  # mu_hat <- matrix(1, ncol = 1, nrow = dim(Delta_phi)[1]) %*% apply(Delta_phi, 2, mean)
  Level_phi <- X[-N,] # removing the last observation to match the dimension with Delta_Y
  R1 <- lm(Level_phi ~ 1)$residuals
  M <- t(R1) %*% R1
  Q <- R1 %*% ginv(M) %*% t(R1)
  decomp <- Re(eigen(t(R0)%*% Q %*% (R0))$values)
  S <- sum((R0 - Q %*% R0)^2)/(N*p-p*l)
  mu <- 2*S*(p+l)

  if(figure){
    plot(decomp, type = "h", xlab = "r", ylab = "")
    abline(h = mu, lty = 2)
  }

  list(r = length(which(decomp >= mu)), lambdas = decomp, mu = mu)
}

rank.buneaAggregated <- function(Z0, Z1, figure = FALSE){
  p <- ncol(Z0)
  N <- nrow(Z1)
  l <- rankMatrix(Z1)[1]

  R0 <- lm(Z0 ~ 1)$residuals
  R1 <- lm(Z1 ~ 1)$residuals
  Minv <- ginv(t(R1) %*% R1)
  # Q <- R1 %*% ginv(M) %*% t(R1)
  # decomp <- Re(eigen(t(R0)%*% Q %*% (R0))$values)
  decomp <- Re(eigen(t(R0)%*% R1 %*% Minv %*% t(R1) %*% (R0))$values)
  # S <- sum((R0 - Q %*% R0)^2)/(N*p-p*l)
  # S <- sum(((diag(x=1, nrow = N) - R1 %*% ginv(M) %*% t(R1) %*% R0)^2)/(N*p-p*l)
  # Frobenius norm is calculated via a trace and a suitable rearrangements - it
  # reduces the dimension of the matrix and hence the computer memory
  S <- sum(diag(t(R0) %*% R0 - t(R0) %*% R1 %*% Minv %*% t(R1) %*% R0 -
                  t(R0) %*% R1 %*% Minv %*% t(R1) %*% R0 +
                  t(R0) %*% R1 %*% Minv %*% t(R1) %*% R1 %*% Minv %*% t(R1) %*% R0))/(N*p-p*l)
  mu <- 2*S*(p+l)
  r.est <- length(which(decomp >= mu))

  if(figure){
    plot(decomp, type = "h", xlab = "r", ylab = "")
    abline(h = mu, lty = 2)
    abline(v = r.est, lty = 2, col = "red")
  }

  list(r = r.est, lambdas = decomp, mu = mu)
}

rank.ic <- function(Z0, Z1, crit = "AIC", dt = 1){
  p <- ncol(Z0)
  N <- nrow(Z1)
  l <- rankMatrix(Z1)[1]

  # ML procedure to get the eigenvalues, r is set arbitrarily to 1
  fitvecm <- vecmAggregated(Z0, Z1, 1, diag(rep(1, p)), diag(rep(1, p)), dt)
  lambdas <- fitvecm[2*1+p+3,]

  if(crit == "AIC"){
    ic <- N*cumsum(log(1-lambdas)) + 2*(1:p)*(2*p-(1:p))+p*(p+1) + 2*p
  } else if(crit == "BIC"){
    ic <- N*cumsum(log(1-lambdas)) + log(N)*((1:p)*(2*p-(1:p))+p*(p+1)/2 + p)
  } else {
    ic <- N*cumsum(log(1-lambdas)) + 2*log(log(N))*((1:p)*(2*p-(1:p))+p*(p+1)/2 + p)
  }

  list(r = which.min(ic), values = ic)
}

rank.ic.alt <- function(Z0, Z1, crit = "AIC", dt = 1){
  p <- ncol(Z0)
  N <- nrow(Z1)
  l <- rankMatrix(Z1)[1]

  # Procedure to get the eigenvalues etc.
  P.Z1 <- Z1 %*% solve(t(Z1) %*% Z1) %*% t(Z1)
  hatPi <- solve(t(Z1) %*% Z1) %*% t(Z1) %*% Z0
  Se <- t(Z0) %*% (diag(rep(1,N)) - P.Z1) %*% Z0
  Sh <- t(hatPi) %*% solve(t(Z1) %*% Z1) %*% hatPi

  lambdas <- c(0, eigen(Sh %*% solve(Se))$values)
  # Zero has been added to cover the case of r=0

  if(crit == "AIC"){
    ic <- N*cumsum(log(1+lambdas)) - 2*(p-(0:(p)))^2
  } else if(crit == "BIC"){
    ic <- N*cumsum(log(1-lambdas)) + log(N)*((1:p)*(2*p-(1:p))+p*(p+1)/2)
  } else {
    ic <- N*cumsum(log(1-lambdas)) + 2*log(log(N))*((1:p)*(2*p-(1:p))+p*(p+1)/2)
  }

  list(r = which.min(ic), values = ic)
}
