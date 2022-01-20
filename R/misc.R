# Miscellaneous tools for analyzing data

# Hilbert transformation using Fourier transforms
    hilbert <- function(inVec,n=length(inVec)){
      #if(!exists("n")) n = length(inVec)
      tmp = fft(inVec[1:n])
      hVec = rep(0,n)
      hVec[c(1,n/2+1)] = c(1,1)
      hVec[2:(n/2)] = 2
      outVec = fft(tmp*hVec,inverse=T)/n

      return(outVec[1:n])
    }

# Unwrap phases function
    unwrapPhase <- function(Z){
      X = as.matrix(Z[,seq(2,ncol(Z),2)]) # Extract x-values
      Y = as.matrix(Z[,seq(3,ncol(Z),2)]) # Extract y-values

      p = ncol(X) # no. of oscillators
      N = nrow(Z) # no. of observations

      phi = matrix(0,nr=N,nc=p)
      for(i in 1:p){
        phi[,i] = signal::unwrap((atan2(Y[,i],X[,i])+2*pi)%%(2*pi))
      }
      phi = data.frame(phi)
      names(phi) = paste0("phi",1:p)
      invisible(phi)
    }

# # Mean square error function
#     mse <- function(x,y){
#       round(mean((x-y)^2),5)
#     }
#

# Mean phase coherence measure (bilateral)
meanPhaseCoherence <- function(phi1,phi2){
  mpc = sqrt(Re(mean(exp(1i*(phi1-phi2))))^2+Im(mean(exp(1i*(phi1-phi2))))^2)
  return(mpc)
}

# # Mean phase coherences
#     mpc <- function(phi,verbose=FALSE){
#       p = ncol(phi)
#       tmp = combn(p,2)
#       Rout = numeric(choose(p,2))
#       for(i in 1:choose(p,2)){
#         Rout[i] = cods::meanPhaseCoherence(phi[,tmp[1,i]],phi[,tmp[2,i]])
#         names(Rout)[i] = paste0("R",tmp[1,i],tmp[2,i])
#       }
#       if(verbose){print(data.frame(mpc=Rout))}
#       invisible(data.frame(mpc=Rout))
#     }
#
# # Check I(1) conditions for the system:
#     checkI1conditions <- function(alpha, beta, verbose=FALSE){
#       alpha
#       p = nrow(alpha)
#       r = ncol(alpha)
#       alphaPerp = qr.Q(qr(alpha), complete = TRUE)[,(r+1):p]
#       betaPerp = qr.Q(qr(beta), complete = TRUE)[,(r+1):p]
#       check1 = abs(sum(t(alphaPerp)%*%alpha)) < 1e-10
#       check2 = abs(sum(t(betaPerp)%*%beta)) < 1e-10
#       check3 = det(t(alphaPerp)%*%betaPerp) != 0
#       if(verbose){
#         print(check1)
#         print(check2)
#         print(check3)
#       } else {
#         if(!check1|!check2|!check3) "I(1) conditions not satisfied"
#       }
#     }

