# TWO sample case simulations
library(otinference)

# Creating a normalized (devided by sigmahat) sample of size K 

  # K = size of the sample 
  # p,q are the given marginals 
  # cost is the cost matrix 
  # n,m are the number of random variables determing the empirical marginals, respectively
  # lambda is the regularization parameter

TwoSinkhornSample <- function(B,p,q,cost,n,m,lambda){
  valueStd <- rep(0,B)
  valueNonStd <- rep(0,B)
  SinkhornTrue <- SinkhornDistance(as.matrix(p),as.matrix(q),as.matrix(cost),lambda)$Distance
  for(i in 1:B){
    phat <- (1/n)*rmultinom(1,n,p)
    qhat <- (1/m)*rmultinom(1,m,q)
    EmpSinkhorn <- SinkhornDistance(as.matrix(phat),as.matrix(qhat),as.matrix(cost),lambda)$Distance
    valueNonStd[i] <- sqrt(m*n/(m+n))*(EmpSinkhorn-SinkhornTrue)
    sigmahat <- (1/sqrt(sigmaTwoCholesky(as.matrix(phat),as.matrix(qhat),n=n,m=n,as.matrix(cost),lambda)))
    valueStd[i] <- sqrt(m*n/(m+n))*sigmahat*(EmpSinkhorn-SinkhornTrue)
  }
  return(data.frame(sampleStd=valueStd, sampleNonStd = valueNonStd,cat = rep("Finite Sample", B)))
}

########### 1. Density Simulation ##############

B <- 10000 #Sample size
alpha <- c(1,1) #Concentration parameter for p and q, respectively
nLambda <- 5  #Lambda = median(cost)/nLambda
L <- 2 #Grid size
n <- 1000
m <- 1000

  # define cost matrix and correct regularization parameter
  cost <- as.matrix(dist(as.matrix(expand.grid(seq(0,1,length.out=L),rev(seq(0,1,length.out=L)))),upper=TRUE,diag=TRUE))
  lambda <- median(cost)/nLambda
  # Define ground measure.
  N <- L^2
  p <- t(as.matrix(rdirichlet(1,rep(alpha[1],N))))
  q <- t(as.matrix(rdirichlet(1,rep(alpha[2],N))))
  
  # Create the sample
  Sample <- TwoSinkhornSample(B,p,q,cost,n,m,lambda)$sampleStd

  #Visualization of the results in comparison to a standard normal density
  ggplot(Sample, aes(x=sample,fill=cat)) + 
    geom_density(alpha=0.1,linetype = "dashed",fill="blue") +
    geom_area(stat = "function", fun = dnorm, fill = "green",alpha=0.1,color="black") +
   # theme_bw(base_size = 24) +
    scale_x_continuous(name = "") +
    scale_y_continuous(name = "Density") +
    guides(fill = FALSE) +
    ggtitle(paste("Sample size n = ", B , " for grid size ",L,"x",L,  sep = ""))
  
########### 2. Kolmogorov-Smirnov Distance Simulation ##############
  
  set.seed(1)
  
  TwoDistance <- data.frame(Lambda = lambda, L = numeric(0), n = numeric(0), m = numeric(0), distance = numeric(0))
  B <- 20000 #Sample size
  alpha <- c(1,1) #Concentration parameter for p and q, respectively
  nLambda <- c(0.01,0.1,1,10) #Lambda = median(cost)/nLambda
  LVals <- c(3,5)
  nVals <- c(50, 100, 150, 500)
  Rep <- 5 # how many measures should be used to average
  
  time <- proc.time()
  for(k in 1:Rep){
    print(paste("Repetition number = ",k, sep=""))
    for(L in LVals){
      print(paste("grid of size = ", L, "x", L, sep=""))
      #Define cost matrix
      cost <- as.matrix(dist(as.matrix(expand.grid(seq(0,1,length.out=L),rev(seq(0,1,length.out=L)))),upper=TRUE,diag=TRUE))
      N <- L^2
      # Define ground measure.
      p <- t(as.matrix(rdirichlet(1,rep(alpha[1],N))))
      q <- t(as.matrix(rdirichlet(1,rep(alpha[2],N))))
      for(lambda in nLambda){
        for(n in nVals){
          print(paste("n = ", n, sep=""))
          print(paste("Lambda = ",signif(median(cost)/lambda,2), sep = ""))
          sample <- TwoSinkhornSample(B,p,q,cost,n,n,median(cost)/lambda)$sampleStd
          Distance <- ks.test(sample,"pnorm")$statistic
          TwoDistance <- rbind(TwoDistance, 
                                      data.frame(Lambda = lambda, L = L, n = n, m = n, distance = Distance))
        }
      }
    }
  }
  print(proc.time()-time)
  
  #Visualization of the results in comparison to a standard normal density
  ggplot(TwoDistance, aes(x=n, y=distance, group=L, color=factor(L))) + 
    stat_summary(aes(color=factor(L)),
                 fun.y=mean, geom = "line",size=2) + 
    stat_summary(aes(colour = factor(L), shape = factor(L)), 
                 fun.y = mean, geom = "point", size = 6) +
    scale_y_log10(limits = c(0.01, max(TwoDistance$distance))) +
    labs(y = "Kolmogorov-Smirnov\n Statistic", x = "Sample size", 
         colour = "L", shape = "L", linetype = "L") +
    scale_x_log10(breaks = unique(TwoDistance$n)) +
    theme_bw(base_size = 2*18) + 
    theme(legend.key.height=unit(0.5,"inch"))+
    facet_wrap(~Lambda)
  
  
########### 3. Comparison to Wasserstein Inference ##############
 
  
  set.seed(1)
  
  TwoDistance <- data.frame(Lambda = numeric(0), L = numeric(0), n = numeric(0), m = numeric(0), distance = numeric(0), comparison = numeric(0))
  B <- 20000 #Sample size
  alpha <- c(1,1) #Concentration parameter for p and q, respectively
  nLambda <- c(0.01,0.1,1,10) #Lambda = median(cost)/nLambda
  LVals <- c(3,5,10)
  nVals <- c(50, 100, 1000, 5000)
  Rep <- 5 # how many measures should be used to average
  
  time <- proc.time()
  for(k in 1:Rep){
    print(paste("Repetition number = ",k, sep=""))
    for(L in LVals){
      print(paste("grid of size = ", L, "x", L, sep=""))
      #Define cost matrix
      cost <- as.matrix(dist(as.matrix(expand.grid(seq(0,1,length.out=L),rev(seq(0,1,length.out=L)))),upper=TRUE,diag=TRUE))
      N <- L^2
      # Define ground measure.
      p <- t(as.matrix(rdirichlet(1,rep(alpha[1],N))))
      q <- t(as.matrix(rdirichlet(1,rep(alpha[2],N))))
      SampleLimit <- limDisAlt(B,c(p),c(q),cost)
      for(lambda in nLambda){
        for(n in nVals){
          print(paste("n = ", n, sep=""))
          print(paste("Lambda = ",signif(median(cost)/lambda,2), sep = ""))
          sample <- TwoSinkhornSample(B,p,q,cost,n,n,median(cost)/lambda)
          DistanceTrue <- ks.test(sample$sampleStd,"pnorm")$statistic
          DistanceComparison <- ks.test(sample$sampleNonStd,SampleLimit)$statistic
          TwoDistance <- rbind(TwoDistance, 
                                   data.frame(Lambda = lambda, L = L, n = n, m = n, distance = DistanceTrue, Comparison = as.factor("No")))
          TwoDistance <- rbind(TwoDistance, 
                                         data.frame(Lambda = lambda, L = L, n = n, m = n, distance = DistanceComparison, Comparison = as.factor("Yes")))
        }
      }
    }
  }
  print(proc.time()-time)
  
  #Visualization of the results for the different regularization parameters 
  ggplot(TwoDistance[(TwoDistance$Comparison=="Yes"),], aes(x=n, y=distance, group=L, color=factor(L))) + 
    stat_summary(aes(color=factor(L)),
                 fun.y=mean, geom = "line",size=2) + 
    stat_summary(aes(colour = factor(L), shape = factor(L)), 
                 fun.y = mean, geom = "point", size = 6) +
    scale_y_log10() +
    labs(y = "Kolmogorov-Smirnov\n Statistic", x = "Sample size", 
         colour = "L", shape = "L", linetype = "L") +
    scale_x_log10(breaks = unique(TwoDistance$n)) +
    theme_bw(base_size = 2*18) + 
    theme(legend.key.height=unit(0.5,"inch")) +
    facet_wrap(~Lambda)
  
  