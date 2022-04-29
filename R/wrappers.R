# Wrapper functions for simulating data,
# performing bootstrap test and inference (restrictions on matrices)
# Output here is more user friendly and usable.


# Wrapper function for the Johansen cpp function.
johansen <- function(X = NULL, r = 1, A = NULL, H = NULL, dt = 0.1, normalize = TRUE){
  X = as.matrix(X)
  N = nrow(X)-1
  p = ncol(X)

  if(is.null(A)){
    A = as.matrix(diag(p)) # produces a unit matrix
  }
  if(is.null(H)){
    H = as.matrix(diag(p))
  }
  df = r*(p-ncol(A))+r*(p-ncol(H))

  out = johansenCpp(X, r, as.matrix(A), as.matrix(H), dt, normalize)

  if(r > 0){
    a.hat = matrix(out[1:r,],nr=p,byrow = TRUE)
    b.hat = matrix(out[(r+1):(2*r),],nr=p,byrow = TRUE)
    rownames(a.hat) = paste0("x",1:p)
    colnames(a.hat) = paste0("r",1:r)
    rownames(b.hat) = paste0("x",1:p)
    colnames(b.hat) = paste0("r",1:r)
  }
  P.hat = out[2*r+1,]
  O.hat = out[(2*r+2):(2*r+1+p),]
  test  = out[(2*r+2+p),]
  eigs  = out[(2*r+2+p+1),]
  res   = out[(2*r+2+p+2):(2*r+2+p+N+1),]

  names(P.hat)    = paste0("x",1:p)
  rownames(O.hat) = paste0("x",1:p)
  colnames(O.hat) = paste0("x",1:p)
  names(test)     = paste0("r=",0:(p-1))
  names(eigs)     = paste0("l",1:p)
  colnames(res)   = paste0("x",1:p)

  if(r==0){
    a.hat = b.hat = rep(0, p)
  }
  return(list(N=N, p=p, r=r ,alpha = a.hat, beta = b.hat, Psi=P.hat, Omega = O.hat,
              test=test, lambda=eigs, A=A, H=H, df=df, dt=dt, res=res, data=X))
}

johansenAggregated <- function(Y = NULL, Z = NULL, r = 1, A = NULL, H = NULL,
                               dt = 0.1, intercept = TRUE, normalize = TRUE){
  N = nrow(Y)
  p = ncol(Y)

  if(is.null(A)){
    A = as.matrix(diag(p)) # produces a unit matrix
  }
  if(is.null(H)){
    H = as.matrix(diag(p))
  }
  df = r*(p-ncol(A))+r*(p-ncol(H))

  out = johansenCppAggregated(Y, Z, r, as.matrix(A), as.matrix(H), dt, intercept, normalize)

  if(r > 0){
    a.hat = matrix(out[1:r,],nr=p,byrow = TRUE)
    b.hat = matrix(out[(r+1):(2*r),],nr=p,byrow = TRUE)
    rownames(a.hat) = paste0("x",1:p)
    colnames(a.hat) = paste0("r",1:r)
    rownames(b.hat) = paste0("x",1:p)
    colnames(b.hat) = paste0("r",1:r)
  }
  P.hat = out[2*r+1,]
  O.hat = out[(2*r+2):(2*r+1+p),]
  test  = out[(2*r+2+p),]
  eigs  = out[(2*r+2+p+1),]
  res   = out[(2*r+2+p+2):(2*r+2+p+N+1),]

  meanY <- apply(Y, 2, mean)

  res0  = Y - matrix(meanY, nrow = N, ncol = p, byrow = T)
  R2    = 1 - sum(res^2)/sum(res0^2)

  names(P.hat)    = paste0("x",1:p)
  rownames(O.hat) = paste0("x",1:p)
  colnames(O.hat) = paste0("x",1:p)
  names(test)     = paste0("r=",0:(p-1))
  names(eigs)     = paste0("l",1:p)
  colnames(res)   = paste0("x",1:p)

  if(r==0){
    a.hat = b.hat = rep(0, p)
  }
  return(list(N = N, p = p, r = r, alpha = a.hat, beta = b.hat, Psi = P.hat, Omega = O.hat,
              test = test, lambda = eigs, A = A, H = H, df = df, dt = dt, r2 = R2))
}

# Wrapper function for the bootstrap cpp function.
bootstrap <- function(X = NULL, B = 1000, dt = 0.1){

  X = as.matrix(X)
  N = nrow(X)-1
  p = ncol(X)

  r = 1
  tmp = johansenCpp(X, r, diag(p), diag(p), dt)
  out = list()
  #
  out$test = tmp[(2*r+2+p),]
  out$boot = bootstrapCpp(X, B = B, dt = dt)
  return(out)
}

bootstrapAggregated <- function(Y = NULL, Z = NULL, B = 1000, dt = 0.1,
                                n = nrow(Y), n.epoch = 1){

  N = nrow(Y)
  p = ncol(Y)

  r = 1
  tmp = johansenCppAggregated(Y, Z, r, diag(p), diag(p), dt)
  out = list()
  #
  out$test = tmp[(2*r+2+p),]
  out$boot = bootstrapCppAggregated(Y, Z, B = B, dt = dt, n = n, n_epoch = n.epoch)
  return(out)
}
