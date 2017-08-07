# ONE sample case simulations
library(otinference)

# Creating a normalized (devided by sigmahat) sample of size K 

# K = size of the sample 
# p,q are the given marginals 
# cost is the cost matrix 
# n is the number of random variables determing the empirical marginal of p
# lambda is the regularization parameter

OneSinkhornSample <- function(B,p,q,cost,n,lambda){
  valueStd <- rep(0,B)
  valueNonStd <- rep(0,B)
  SinkhornTrue <- SinkhornDistance(as.matrix(p),as.matrix(q),as.matrix(cost),lambda)$Distance
  for(i in 1:B){
    phat <- (1/n)*rmultinom(1,n,p)
    EmpSinkhorn <- SinkhornDistance(as.matrix(phat),as.matrix(q),as.matrix(cost),lambda)$Distance
    valueNonStd[i] <- sqrt(n)*(EmpSinkhorn-SinkhornTrue)
    sigmahat <- (1/sqrt(sigmaOneCholesky(as.matrix(phat),as.matrix(q),as.matrix(cost),lambda)))
    valueStd[i] <- sqrt(n)*sigmahat*(EmpSinkhorn-SinkhornTrue)
   # print(i)
  }
  return(data.frame(sampleStd=valueStd, sampleNonStd = valueNonStd, cat = rep("Finite Sample", B)))
}

########### 1. Density Simulation ##############

B <- 20000 #Sample size
alpha <- c(1,1) #Concentration parameter for p and q, respectively
nLambda <- 1  #Lambda = median(cost)/nLambda
L <- 10 #Grid size
n <- 1000

# define cost matrix and correct regularization parameter
cost <- as.matrix(dist(as.matrix(expand.grid(seq(0,1,length.out=L),rev(seq(0,1,length.out=L)))),upper=TRUE,diag=TRUE))
lambda <- median(cost)/nLambda
# Define ground measure.
N <- L^2
p <- t(as.matrix(rdirichlet(1,rep(alpha[1],N))))
q <- t(as.matrix(rdirichlet(1,rep(alpha[2],N))))

# Create the sample
Sample <- OneSinkhornSample(B,p,q,cost,n,lambda)$sampleStd

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

OneDistance <- data.frame(Lambda = numeric(0), L = numeric(0), n = numeric(0), distance = numeric(0), Comparison = numeric(0))

B <- 500 #Sample size
alpha <- c(1,1) #Concentration parameter for p and q, respectively
nLambda <- c(0.01,0.1,1,10) #Lambda = median(cost)/nLambda
LVals <- c(3,5)
nVals <- c(50,100,1000,5000)
Rep <- 2 # how many measures should be used to average

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
    q <- p
    SampleLimit <- limDisNull(B,c(p),cost)
    for(lambda in nLambda){
      for(n in nVals){
        print(paste("n = ", n, sep=""))
        print(paste("Lambda = ",signif(median(cost)/lambda,2), sep = ""))
        sample <- OneSinkhornSample(B,p,q,cost,n,median(cost)/lambda)
        DistanceTrue <- ks.test(sample$sampleStd,"pnorm")$statistic
        DistanceComparison <- ks.test(sample$sampleNonStd,SampleLimit)$statistic
        OneDistance <- rbind(OneDistance, 
                             data.frame(Lambda = lambda, L = L, n = n, distance = DistanceTrue, Comparison = as.factor("No")))
        OneDistance <- rbind(OneDistance, 
                                       data.frame(Lambda = lambda, L = L, n = n, distance = DistanceComparison, Comparison = as.factor("Yes")))
      }
    }
  }
}
print(proc.time()-time)

#Visualization of the results for the different regularization parameters 
ggplot(OneDistance[(OneDistance$Comparison=="No"),], aes(x=n, y=distance, group=L, color=factor(L))) + 
  stat_summary(aes(color=factor(L)),
               fun.y=mean, geom = "line",size=2) + 
  stat_summary(aes(colour = factor(L), shape = factor(L)), 
               fun.y = mean, geom = "point", size = 6) +
  scale_y_log10(limits = c(0.01, max(OneDistance$distance))) +
  labs(y = "Kolmogorov-Smirnov\n Statistic", x = "Sample size", 
       colour = "L", shape = "L", linetype = "L") +
  scale_x_log10(breaks = unique(OneDistance$n)) +
  theme_bw(base_size = 2*18) + 
  theme(legend.key.height=unit(0.5,"inch"))+
  facet_wrap(~Lambda)
