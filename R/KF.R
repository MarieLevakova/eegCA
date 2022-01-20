# MULTIVARIATE REGRESSION CRITERION
# Author: Christophe Giraud, Ecole Polytechnique, FRANCE
# contact: christophe.giraud@polytechnique.edu
# date: 28/10/2010

KF <- function(Y,X,rmax=min(dim(X)[1],dim(X)[2],dim(Y)[2]),K=2,CV=FALSE,V=10,rkX=0,tol=10**(-6)){
	# INPUT
	# Y: m x n matrix
	# X: m x p matrix
	# 
	# OUPUT
	# XA: m x n matrix
	# A: p x n matrix
	# rank: integer
	#
	m<-dim(Y)[1]
	n<-dim(Y)[2]
	p<-dim(X)[2]
	if (rkX==0) rkX <- effectiverank(X,tol)
	rmax <- min(rmax,rkX)
	if (CV) {
		K=(11:30)/10
		penK<-pen(n,m,rkX,rmax,K)
		mv=round(m/V)
		rmaxv=min(rmax,m-mv)
		penCV <- pen(n,m-mv,min(m-mv,rkX),rmaxv,K)
		res <- hatACV(Y,X,rmax,K=K,penK=penK,penCV=penCV,V=V,rkX=rkX)
	}
	else {
		penK<-pen(n,m,rkX,rmax,K)
		res <- hatA(Y,X,rmax,K=K,penK=penK,rkX=rkX)
	}
	dimnames(res$XA) <- list(paste("K=",K),NULL,NULL)
	dimnames(res$A)  <- list(paste("K=",K),NULL,NULL)
	return(list(XA=drop(res$XA),A=drop(res$A),rank=drop(res$rank)))
}

hatACV <- function(Y,X,rmax,K,penK,penCV,V,rkX){ 
	m<-dim(Y)[1]
	n<-dim(Y)[2]
	p<-dim(X)[2]
	mv <- round(m/V)
	errorV <-array(0,c(V,length(K)))
	for (iv in 1:V){
		It=(((iv-1)*mv+1):(iv*mv))
		Ia=(1:m)[-It]
		Xa <- X[Ia,]
		Ya <- Y[Ia,]
		chapAK <- hatA(Ya,Xa,min(rmax,m-mv),K,penCV,min(rkX,m-mv))$A
		errorV[iv,]<-errorCV(Y[It,],X[It,],chapAK,K)
		}
	error <- apply(errorV,2,mean)
	iK <- which.min(error)
	KCV <- K[iK]
	res <- hatA(Y,X,rmax,K=c(KCV,1),penK=penK[c(iK,1),],rkX)
	XchapACV <- res$XA[1,,]
	rankCV <- res$rank[1]
	return(list(XA=XchapACV,A=res$A[1,,],rank=rankCV))	
}


hatA <- function(Y,X,rmax,K,penK,rkX){
	m<-dim(Y)[1]
	n<-dim(Y)[2]
	p<-dim(X)[2]
	k<-length(K)
	rmax<-min(c(rmax,p,n,m,rkX))
	chapXAr <- hatAr(Y,X,rmax=rmax,rkX=rkX)$XA
	crit <- array(0,c(rmax,k))
	for (r in 1:rmax) {
		crit[r,] <-	sum((Y-chapXAr[r,,])**2)*(1+penK[,r])
	} 
	rchap<- apply(crit,2,which.min)
	chapXA <- array(0,c(k,m,n))
	chapA <- array(0,c(k,p,n))
	NormY <- sum(Y**2)
	for (ik in 1:k) {
		chapXA[ik,,]<-chapXAr[rchap[ik],,]
		if (NormY < crit[rchap[ik],ik]) {
			chapXA[ik,,]<-chapXA[ik,,]*0
			rchap[ik] <- 0
		}
		chapA[ik,,]<-inverseM(t(X)%*%X,rkX=rkX)%*%t(X)%*%chapXA[ik,,]
		}
	return(list(XA=chapXA,rank=rchap,A=chapA))
}


hatAr <- function(Y,X,rmax,rkX,computeA=FALSE){
	# Y: mxn
	# X: mxp
	#rmax: integer
	m<-dim(Y)[1]
	n<-dim(Y)[2]
	p<-dim(X)[2]
	rmax<-min(rmax,min(n,m),rkX)
	Ar<-array(0,c(rmax,p,n))
	XAr<-array(0,c(rmax,m,n))
	PY <- X%*%inverseM(t(X)%*%X,rkX=rkX)%*%t(X)%*%Y
	svdPY<-svd(PY,nu=rmax,nv=rmax)
	for (r in 1:rmax){
		if (r==1) XAr[r,,]<-svdPY$u[,1:r]%*%t(svdPY$v[,1:r])*svdPY$d[1]
		else 	XAr[r,,]<-svdPY$u[,1:r]%*%diag(svdPY$d[1:r])%*%t(svdPY$v[,1:r])
		if (computeA) 	Ar[r,,]<-inverseM(t(X)%*%X,rkX=rkX)%*%t(X)%*%XAr[r,,]
	}
	return(list(A=Ar,XA=XAr))
}

effectiverank <- function(X,tol){
	singular <- svd(X,nu=0,nv=0)$d
	rkX <- sum(sign(singular[singular>tol]))
	return(rkX)
}


pen <- function(n,m,q,rmax,K,Nsim=200){
	# n,m,rmax,Nsim,q: integer
	# K: vector
	penK <- array(0,c(length(K),rmax))
	if (q*n<1000) Sr2 <- (SMC(q,n)[1:rmax])**2
	else Sr2 <- (Sapprox(q,n)[1:rmax])**2	
	for (iK in 1:length(K)){
		penK[iK,]<- K[iK]*Sr2/pmax(n*m-1-K[iK]*Sr2,10**(-10))
	}
	dimnames(penK) <- list(paste("K=",K),NULL)		
	return(penK) # length(K) x rmax
}

SMC <- function(q,n,Nsim=200){
	Sk <- array(0,c(Nsim,min(q,n)))	
		for (is in 1:Nsim) {
			s <- svd(matrix(rnorm(q*n),nrow=q,ncol=n),nu=0,nv=0)$d
			Sk[is,]<-sqrt(cumsum(s**2))
	}
	return(apply(Sk,2,mean))		
}

Sapprox <-function(q,n,eps=10**(-9)){
	beta <- min(n,q)/max(n,q)
	alpha <- (1:min(n,q))/min(n,q)
	s<-rep(0,min(n,q))
	f <- function(x){return(sqrt((x-(1-sqrt(beta))^2)*((1+sqrt(beta))^2-x))/(2*pi*beta*x))}
	xf <- function(x){return(sqrt((x-(1-sqrt(beta))^2)*((1+sqrt(beta))^2-x))/(2*pi*beta))}
	for (a in 1:(length(alpha)-1)){
		m <- (1-sqrt(beta))^2
		M <- (1+sqrt(beta))^2
		while ((M-m)>eps) {
			if (integrate(f,(m+M)/2,(1+sqrt(beta))^2,rel.tol=eps/100)$value<alpha[a]) M <- (m+M)/2
			else m <- (m+M)/2
		}
		s[a] <- integrate(xf,(m+M)/2,(1+sqrt(beta))^2,rel.tol=eps/100)$value
	}
	s[length(alpha)] <- 1
	return(sqrt(s*n*q))
} 


errorCV<- function(Yt,Xt,chapAl,l){
	ll=length(l)
	error <- rep(0,ll)  
	for (il in 1:ll)	error[il] <- sum((Yt-Xt%*%chapAl[il,,])**2)
	return(error)
	}


inverseM <- function(M,rkX,tol=10**(-7)){
	svdM <- svd(M,nu=rkX,nv=rkX)
	singul <- svdM$d[1:rkX]	
	Itol <- (1:length(singul))[singul>tol]
	invM <- array(0,dim(M))
	if (length(Itol)>0)	invM <- svdM$v[,Itol]%*%diag(1/svdM$d[Itol],nrow=length(Itol))%*%t(svdM$u[,Itol])
	return(invM)
}

generateData <- function(m,n,p,rho,b,r){
	Y <- array(0,c(m,n))
	Sigma<- array(0,c(p,p))
	for (j in 1:p) Sigma[j,]<-rho**(abs(j-1:p))
	e<-eigen(Sigma,symmetric=TRUE)	
	B0<-matrix(rnorm(p*r),ncol=r)
	B1<-matrix(rnorm(r*n),ncol=n)
	A<-b*B0%*%B1 # A est pxn
	X <- array(0,c(m,p))
	for (i in 1:m)  {
        X[i,]<- e$vectors%*%diag(sqrt(e$values))%*%e$vectors%*%rnorm(p)
		Y[i,] <- t(A)%*%X[i,]+rnorm(n)
	}
	return(list(Y=Y,A=A,X=X))
}


