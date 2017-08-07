#Computation for sigmahat in the one and two sample case, respectively. The computation is performed via sparse matrices and cholesky decomposition.

library(Matrix) #for sparse matrix and the command bdiag


########################################################## OneSampleCase ########################################

sigmaOneCholesky <- function(p,q,cost,lambda){
  
  #First, we select the support of p and q
  cost <- as.matrix(cost[p>0,q>0])
  p <- as.matrix(p[p>0])
  q <- as.matrix(q[q>0])
  
  #1.Covariance structure
  p_temp <- drop(p)
  varp <- p_temp*(1-p_temp)
  Covariance <- -outer(p_temp, p_temp)
  diag(Covariance) <- varp
  
  #2.Delta Pi, i.e. we require the coefficient matrix A and the inverse of the second derivative of the Lagrangian evaluated at the optimal solution denoted by LagraInv
  d <- length(p)
  h <- length(q)
  #Coefficient matrix as sparse matrix generated via the kronecker product
  Coeff1 <- kronecker(Matrix(diag(length(p)),sparse = TRUE),t(as.matrix(rep(1,length(q)))))
  Coeff2 <- kronecker(as.matrix(t(rep(1,length(p)))),Matrix(diag(length(q)),sparse=TRUE))
  Coeff2 <- Coeff2[-dim(Coeff2)[1],]
  Coeff <- Matrix(rbind(Coeff1,Coeff2),sparse=TRUE)
  
  #Lagrangian Inverse
  Transport <- SinkhornDistance(p,q,cost,lambda)$Transportplan  #calculating optimal solution with Rcpp Armadillo
  LagraInv <- Diagonal(length(t(Transport)),t(Transport))
  #LagraInv <- Matrix(0,dim(Transport)[1]*dim(Transport)[2],dim(Transport)[1]*dim(Transport)[2])
  #diag(LagraInv) <- t(Transport)
  
  #Cost vector
  Cost <- matrix(cost,1,d*h)
  
  
  #Delta matrix, i.e., the derivative of the parametrization
    # DA^T 
  A <- t(Coeff)
  A@x <- A@x*matrix(t(Transport),ncol=1)[A@i+1]
    # A(DA^T)
  inverse <- Coeff%*%A  
    # (ADA^T)^(-1)
  inverse <- as.matrix(chol2inv(chol(inverse)))
    # DA^T((ADA^T)^(-1)) (note that A is now equal to DA^T by the first step)
  Delta <- as.matrix(A%*%inverse[,1:d])
  
  #Sigma calculation
  sigma <- as.double(Cost%*%Delta%*%Covariance%*%t(Delta)%*%t(Cost))
  
  return(sigma)
}

########################################################## TwoSampleCase ########################################

sigmaTwoCholesky <- function(p,q,n,m,cost,lambda){
  
  #First, we select the support of p and q
  cost <- as.matrix(cost[p>0,q>0])
  p <- as.matrix(p[p>0])
  q <- as.matrix(q[q>0])
  
  #1.Covariance structure
  p_temp <- drop(p)
  varp <- p_temp*(1-p_temp)
  sigmaP <- -outer(p_temp, p_temp)
  diag(sigmaP) <- varp
  sigmaP <- (m/(n+m))*sigmaP
  q_temp <- drop(q)
  q_temp <- q_temp[-length(q_temp)]
  varq <- q_temp*(1-q_temp)
  sigmaQ <- -outer(q_temp,q_temp)
  diag(sigmaQ) <- varq
  sigmaQ <- (1-m/(n+m))*sigmaQ
  Covariance <- bdiag(sigmaP,sigmaQ)
  
  #2.Delta Pi, i.e. we require the coefficient matrix A and the inverse of the second derivative of the Lagrangian evaluated at the optimal solution denoted by LagraInv
  d <- length(p)
  h <- length(q)
  
  #Coefficient matrix as sparse matrix generated via the kronecker product
  Coeff1 <- kronecker(Matrix(diag(length(p)),sparse = TRUE),t(as.matrix(rep(1,length(q)))))
  Coeff2 <- kronecker(as.matrix(t(rep(1,length(p)))),Matrix(diag(length(q)),sparse=TRUE))
  Coeff2 <- Coeff2[-dim(Coeff2)[1],]
  Coeff <- Matrix(rbind(Coeff1,Coeff2),sparse=TRUE)
  
  #Lagrangian Inverse
  Transport <- SinkhornDistance(p,q,cost,lambda)$Transportplan  #calculating optimal solution with Rcpp Armadillo
  LagraInv <- Diagonal(length(t(Transport)),t(Transport))
  #LagraInv <- Matrix(0,dim(Transport)[1]*dim(Transport)[2],dim(Transport)[1]*dim(Transport)[2])
  #diag(LagraInv) <- t(Transport)
  
  #Cost vector
  Cost <- matrix(cost,1,d*h)
  
  #Delta matrix, i.e., the derivative of the parametrization
    # DA^T
  A <- t(Coeff)
  A@x <- A@x*matrix(t(Transport),ncol=1)[A@i+1]
    # A(DA^T)
  inverse <- Coeff%*%A  
    # (ADA^T)^(-1)
  inverse <- as.matrix(chol2inv(chol(inverse)))
    # DA^T((ADA^T)^(-1)) (note that A is now equal to DA^T by the first step)
  Delta <- as.matrix(A%*%inverse)
  
  #Sigma calculation
  sigma <- as.double(Cost%*%Delta%*%Covariance%*%t(Delta)%*%t(Cost))
  
  return(sigma)
}