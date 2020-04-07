#########################
### CODE FROM MATSUI ####
#########################

### Modified for L1 project: 
### only needed functions kept
### fixed implementation of CalcNormType for type=3 (Gertheiss norm)

### See paper: Matsui - Variable and boundary selection for functional data via multiclass logistic regression modeling (2014)


# calculate generalized inverse
ginv = function (X, tol = sqrt(.Machine$double.eps)) {
  if (length(dim(X)) > 2L || !(is.numeric(X) || is.complex(X))) 
    stop("'X' must be a numeric or complex matrix")
  if (!is.matrix(X)) 
    X <- as.matrix(X)
  Xsvd <- svd(X)
  if (is.complex(X)) 
    Xsvd$u <- Conj(Xsvd$u)
  Positive <- Xsvd$d > max(tol * Xsvd$d[1L], 0)
  if (all(Positive)) 
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive)) 
    array(0, dim(X)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * 
                                               t(Xsvd$u[, Positive, drop = FALSE]))
}

### derivative of L1 penalty
# [in] beta: parameter estimator(p x 1)
# [in] lambda:regularization parameter
# [in] a:tuning parameter (for SCAD, elnet)
# [in] b0:intercept 1/0
# return penalty vector (p x 1)
# date	2012-06-21
L1Pen <- function(beta, lambda, a, b0=T){
	l1pen<-1:length(beta)
	for(j in 1:length(beta))
	{
		l1pen[j]<- a * lambda + (1-a) * lambda * beta[j] #Friedman et al.(2010)
	}
	return(l1pen)
}

## penalty
Pena2 <- function(m){ #return m x m matrix
  D <- matrix(0,nrow=m-2,ncol=m)
  for(i in 1:(m-2)){
    D[i,i] <- 1
    D[i,i+1] <-  -2
    D[i,i+2] <- 1
  }
  DD <- matrix(0,nr=m,nc=m)
  return(t(D)%*%D)
}

### select norm type
# [in] x:vector
# [in] type:type of norm: 2=Aguilera, 3=Gertheiss, else=Euclidean
# [in] lam2 :Gerthesis-type regularization parameter
# date  2014-01-07
CalcNormType = function(x, type, lam2 = 1,m_Jphi=NULL,m_Jphi2=NULL){
  res = 0
  x = as.matrix(x)
  if(nrow(x)<4) type = 1
  if(type == 2){  #Aguilera
    Mat = Pena2(nrow(x))
    res = t(x) %*% Mat %*% x
    res = sqrt(sum(diag(res)))  #if x is a matrix
    
  } else if (type == 3){  #Gerthesis
    Mat = Pena2(nrow(x))
    #res = t(x) %*% (m_Jphi + lam2 * Mat) %*% x # this is not Gerthesis, but this is used in Matzui simulations
    res = t(x) %*% (m_Jphi + lam2 * m_Jphi2) %*% x
    res = sqrt(sum(diag(res)))  #if x is a matrix
    
  } else { #Euclid
    res = sqrt(sum(x^2))
  }
  return (res)
}

### Hessian matrix
# [in] l:variable 1
# [in] k:variable 2
# [in] X:design matrix(n x p)
# [in] prob:probability vector(n x 1)
# return Hessian matrix (p x p)
# date	2012-07-05
ddl_fnc<-function(l,k,X,prob){
	V = W = 1:nrow(X)
	VB = WB = matrix(0, nr=nrow(X), nc=ncol(X))
	for(i in 1:nrow(X)){
		V[i]<-prob[i,l]*prob[i,k]
		W[i]<-prob[i,l]*(1-prob[i,k])
		VB[i,]<-as.vector(V[i]*X[i,])
		WB[i,]<-as.vector(W[i]*X[i,])
	}
	if(l==k){return(-t(X)%*%WB)}
	if(l!=k){return(t(X)%*%VB)}
}

#############################################################################
# Multiclass functional logistic + group lasso (between variable only) estimated by LQA
# Penalty P1
# Require m_p, m_M
# [in] X: predictors (design matrix)
# [in] y: response
# [in] groups: grouping structure (groups of predictors, corresponding to different features)
# [in] lambda: regularization parameter
# [in] a=1: tuning parameter (elastic net) a=1:lasso
# [in] b0: intercept
# [in] type: norm type
# [in] lam2: regularization parameter in CalcNormType (if type=3)
# return parameter estimator
# date: 2013-01-11
#       2014-01-07
#       2018-05-10
gLQA_mlogistic2_mod <- function(X, y, groups, lambda, type=1, lam2=0, a=1, b0=TRUE, m_Jphi=NULL, m_Jphi2=NULL){
	X<-as.matrix(X)
	y<-as.matrix(y)
	K<-ncol(y)+1 # number of classes
	N<-0
	# mvec = rep(m_M, m_p) # number of predictors corresponding to each functional variable
	mvec = as.numeric(table(groups)) # number of predictors corresponding to each variable
	# p = m_p	# number of functional variables
	p = length(mvec) # number of variables
	M1 = b0 + c(1,cumsum(mvec)[-p]+1)	# first predictor for each group (for each variable)

	beta<-matrix(5,nr=ncol(X),nc=K-1)
	#	nbeta<-matrix(0.0001,nr=ncol(X),nc=K-1) # initial value
	nbeta<-matrix(ginv(t(X)%*%X)%*%t(X)%*%y,nr=ncol(X),nc=K-1) # initial value

	prob<-matrix(0,nr=nrow(X),nc=K-1) # column: probability of belonging to each class
	dl<-rep(0,ncol(X)*(K-1))
	ddl<-matrix(0,nr=ncol(X)*(K-1),nc=ncol(X)*(K-1))

	erase = matrix(1,nr=ncol(X), nc=K-1)
	#	eps = 1e-8 / (2*nrow(X)*L1Pen(0,lambda[1]))	# Zou et al (2007) better

	## repetitive start
	while(max(abs(nbeta-beta))>1e-5 && N < 200){
		beta<-nbeta

		# penalty term calculation
		pena<-matrix(0, nr=sum(mvec)+b0, nc=K-1)
		for(l in 1:(K-1)) {
			for(j in 1:p) { # predictor j
        elbuf = M1[j]:(M1[j]+mvec[j]-1)  # elements corresponding to predictor j
        # buf <- L1Pen(sqrt(sum(beta[elbuf, l]^2)), lambda, a, b0) / (sqrt(sum(beta[elbuf, l]^2)) + 0)
        buf <- L1Pen(CalcNormType(beta[elbuf, l], type, lam2, m_Jphi=m_Jphi, m_Jphi2=m_Jphi2), lambda, a, b0) / (CalcNormType(beta[elbuf, l], type, lam2, m_Jphi=m_Jphi, m_Jphi2=m_Jphi2) + 0) # 2014-01-07 change
        pena[elbuf, l] = rep(buf, mvec[j]) * sqrt(mvec[j])	# Square root of number of group elements required
			}
		}
		#	pena <- ifelse(pena==Inf,0,pena)	# With an estimated amount of 0, there is no penalty 0 or not
		Pena <- diag(as.vector(pena))

		for(i in 1:nrow(X)){	# create prob
			for(l in 1:(K-1))
				prob[i,l]<-exp(X[i,]%*%beta[,l]) / (1+sum(exp(X[i,]%*%beta)))
		}

		# Log-likelihood calculation for iteration judgment
		#loglik<-sum(y * log(prob)) + 
		#	sum((1 - apply(y,1,sum)) * log(1 - apply(prob,1,sum)))
		#if(is.na(loglik)) return(list(beta=nbeta,erase=erase))


		for(l in 1:(K-1)) # Create Log-Likelihood first order derivative
			dl[ncol(X)*(l-1)+(1:ncol(X))]<-as.vector(t(y[,l]-prob[,l])%*%X)
		for(l in 1:(K-1)) # create Log likelihood second order differential
			for(k in 1:(K-1))
				ddl[ncol(X)*(l-1)+(1:ncol(X)),ncol(X)*(k-1)+(1:ncol(X))]<-ddl_fnc(l,k,X,prob)
		# calculate updated value
		erase.no = NULL	  # serial number of excluded beta
		for(l in 1:(K-1)) {
			for(j in 1:p) {
				if(erase[M1[j], l] == 0) {
					erase.no =append(erase.no, (M1[j]:(M1[j]+mvec[j]-1)) + (l-1)*nrow(beta))
				}
			}
		}
		# calculate with excluded parts only
		if(!is.null(erase.no)) {
			nbeta[] = 0 
			nbeta[-erase.no]<-as.vector(beta)[-erase.no] -
					solve(ddl[-erase.no,-erase.no] - 
					nrow(X) * Pena[-erase.no,-erase.no])%*%
					(dl[-erase.no] - 
					nrow(X) * Pena[-erase.no, -erase.no] %*% 
					as.vector(beta)[-erase.no])
		} else {
			nbeta[,]<-as.vector(beta) - 
						solve(ddl - nrow(X) * Pena) %*% 
						(dl - nrow(X) * Pena %*% as.vector(beta))
		}
		for(l in 1:(K-1)) {
			for(j in 1:p) { # predictor j
			  # change condition (divide by number) 2013-01-17
			  # if(sum(nbeta[(M1[j]:(M1[j]+mvec[j]-1)), l]^2)/m_M <1e-6) {
			  if(sum(nbeta[(M1[j]:(M1[j]+mvec[j]-1)), l]^2)/mvec[j] <1e-6) {
					#||beta|| If it is small enough remove that variable
					nbeta[(M1[j]:(M1[j]+mvec[j]-1)), l]<-0
					erase[(M1[j]:(M1[j]+mvec[j]-1)), l] = 0
				}
			}
		}
		#browser()
		N<-N+1
	}
	return(list(beta=nbeta,erase=erase))
}

### BIC for group SCAD, Gauss (Wang et al., 2007)
### need to estimate coefficients in advance (function gLQA_mlogistic2_mod)
# [in] X: predictors (design matrix)
# [in] y: response
# [in] groups: grouping structure (groups of predictors, corresponding to different features)
# [in] lambda: regularization parameter
# [in] est: estimation from gLQA_mlogistic2_mod
# [in] a: tuning parameter
# [in] b0: intercept
# [in] type: norm type
# [in] lam2: regularization parameter in CalcNormType (if type=3)
# return FBIC
# date:2013-01-17  
#      2014-01-07
#      2018-05-10
BIC_multilog2_L1_mod <-function(X, y, groups, lambda, beta, type=1, lam2=0, a=1, b0=TRUE, m_Jphi=NULL, m_Jphi2=NULL){
  X<-as.matrix(X)
  y<-as.matrix(y)
  K<-ncol(y)+1 # number of classes
  N<-0
  # mvec = rep(m_M, m_p) # number of predictors corresponding to each functional variable
  mvec = as.numeric(table(groups)) # number of predictors corresponding to each variable
  # p = m_p	# number of functional variables
  p = length(mvec) # number of variables
  M1 = b0 + c(1,cumsum(mvec)[-p]+1)	# first predictor for each group (for each variable)
  
  prob<-matrix(0,nr=nrow(X),nc=K-1) # column: probability of belonging to each class
  n <- nrow(X)	# sample size
  
	# penalty term calculation
	pena<-matrix(0, nr=sum(mvec)+b0, nc=K-1)
	for(l in 1:(K-1)) {
		for(j in 1:p) { # variable j
		  elbuf = M1[j]:(M1[j]+mvec[j]-1)  # elements corresponding to predictor j
		  # buf <- L1Pen(sqrt(sum(beta[elbuf, l]^2)), lambda, a, b0) / (sqrt(sum(beta[elbuf, l]^2)) + 0)
		  buf <- L1Pen(CalcNormType(beta[elbuf, l], type, lam2, m_Jphi=m_Jphi, m_Jphi2=m_Jphi2), lambda, a, b0) / (CalcNormType(beta[elbuf, l], type, lam2, m_Jphi=m_Jphi, m_Jphi2=m_Jphi2) + 0) # 2014-01-07 change
		  pena[elbuf, l] = rep(buf, mvec[j]) * sqrt(mvec[j])	# Square root of number of group elements required
		}
	}
	#	pena <- ifelse(pena==Inf,0,pena)	# With an estimated amount of 0, there is no penalty 0 or not
	Pena <- diag(as.vector(pena))
  
	for(i in 1:nrow(X)){	# create prob
	  for(l in 1:(K-1))
	    prob[i,l]<-exp(X[i,]%*%beta[,l]) / (1+sum(exp(X[i,]%*%beta)))
	}

	# avoid diverging likelihood
	buf = 1 - apply(prob,1,sum)
	for(i in 1:length(buf)) {
		if (buf[i] == 1) buf[i] = 1 - 1e-7
		else if (buf[i] < 1e-15) buf[i] = 1e-15
	}
  
	# Wmat
	W = matrix(0, nr=nrow(X)*(K-1), nc=nrow(X)*(K-1)) #weight matrix p*(1-p)
	for(l in 1:(K-1)){
	  diag(W[(1:n)+(l-1)*n, (1:n)+(l-1)*n]) = prob[, l] * (1 - prob[, l])
	  if (l==1) next
	  for(ll in (l-1):1){
	    diag(W[(1:n)+(ll-1)*n, (1:n)+(l-1)*n]) = -prob[, l] * prob[, ll]
	    diag(W[(1:n)+(l-1)*n, (1:n)+(ll-1)*n]) = -prob[, l] * prob[, ll]
	  }  
	}
	
  # likelihood
	loglik<-sum(y * log(prob)) + sum((1 - apply(y,1,sum)) * log(buf))
	if(is.na(loglik)) return(Inf)
	
	#calculate active Set
  erase = matrix(1,nr=ncol(X), nc=K-1)
  for(l in 1:(K-1)) {
    for(j in 1:p) { # variable j
  	  # eliminated variables
  	  if (sum(abs(beta[M1[j], l])) == 0)
  	    erase[M1[j]:(M1[j]+mvec[j]-1), l] = 0
  	}
  }
  
	# degrees of freedom
	XX = diag(1, K-1) %x% X
	penamat = diag(as.vector(pena))
	
	buf = t(XX) %*% W %*% XX + n*penamat
	active = which(erase==1) #as vector
	if (length(active)==b0*(K-1)) return(Inf) 
	buf = buf[active, active] #active set
	#hat matrix
	Hatmat = W %*% XX[,active] %*% ginv(buf) %*% t(XX[,active])
	
	#BIC
#	print(c(-2*loglik, log(n)*(sum(beta!=0)-b0)))
#	BIC = -2 * loglik + log(n) * (sum(beta!=0)-b0)
#	print(c(-2*loglik, log(n)*sum(diag(Hatmat))))
	BIC = -2 * loglik + log(n) * sum(diag(Hatmat))
	return (BIC)
}
