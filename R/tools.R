# Tools for the package, tests and calcuations.

# Find the p-values for the Johansen procedure using bootstrapped values
    rankTest <- function(M){
      p = length(M$test)
      pVals = numeric(0)
      for(i in 1:p){
        test.cdf = ecdf(M$boot[,i])
        pVals[i] = 1-test.cdf(M$test[i])
      }
      names(pVals) = paste0("r=",(0:(p-1)))
      if(sum(pVals>0.05)!=0){
        r.est = min(which(pVals>0.05))-1
      } else{
        r.est = p
      }
      names(r.est) = ""
      out = list(r=r.est,pVal = pVals)
      return(out)
    }


# Test restrictions (using LRT) on alpha and beta, input models (from:unrestricted/previous, to:restrictions)
    LRtest <- function(from,to){
      T     = nrow(from$res)+1
      r     = to$r
      df    = to$df
      lp    = from$lambda
      lh    = to$lambda
      test  = T*sum(log(1-lh[1:r])-log(1-lp[1:r]))
      pVal  = pchisq(test,df,lower.tail = FALSE)
      out   = round(data.frame(test,df,pVal),3)
      names(out) = c("Statistic","df","p-value")
      return(out)
    }


# Caculate the standard errors for alpha and mu
    getSE <- function(fit){
      N   = fit$N
      dt  = fit$dt
      Om  = fit$Omega
      phi = fit$data
      r   = fit$r
      p   = nrow(Om)

      # Calculate...
      if(r==0){
        SE = sqrt(diag(Om)*dt)
        alpha = NULL
      } else {
        alpha = fit$alpha
        beta = fit$beta
        #Z = t(beta)%*%t(phi)
        Z = crossprod(beta,t(phi))
        Z = rbind(Z,1)
        #ZZ = solve(Z%*%t(Z))
        ZZ = solve(tcrossprod(Z,Z))
        SE = sqrt(diag(kronecker(Om,ZZ))/dt)
      }

      # Set output with headers
      outMat = matrix(SE,nr=3,byrow=TRUE)
      outVec = c(outMat)
      rownames(outMat) = c("phi_1","phi_2","phi_3")
      if(r > 0){
        colnames(outMat)=c(paste0("alpha_",1:r),"mu")
        names(outVec) = c(paste0("alpha_",as.character(outer(10*(1:p),(1:r),FUN = "+"))),paste0("mu_",1:p))
      } else {
        colnames(outMat)=paste0("mu")
        names(outVec)=paste0("mu_",1:p)
      }
      out = rbind(c(alpha,fit$Psi),outVec)
      tVal = abs(out[1,]/out[2,])
      pVal = 2*pnorm(tVal,lower.tail = FALSE)
      out = rbind(out,pVal)

      colnames(out) = names(outVec)
      rownames(out) = c("Estimate","Std.Err","p-value")


      # Tvals = abs(estPar)/stdErr
      # pVals = 2*pnorm(Tvals,lower.tail = FALSE)
      return(t(out))
    }

