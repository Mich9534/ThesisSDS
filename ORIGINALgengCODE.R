rm(list = ls())

setwd('/home/michele/Desktop/Small_code_R')
source("cogibbs_gn.R")
source("cogibbs_po.R")
# source("fgenera_data.R")
source('Vn_Miller.R')
source('fv_it.R')
source('fgamk_i.R') # here i = n already implemented into fgamk
source('fk_t.R')  
source('ApproximationVn.R')
source('lVnGn.R')
source('vec2mat.R')
source('pr_cc.R')
source('edge_est.R')

library(doSNOW)
require(parallel)
library(gdata)
library(extraDistr)
library(fossil)
library(reshape)
library(gdata)
library(igraph)
library(mcclust.ext)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(coda)
library(dummies)
library(randnet) 
# library(igraphdata)
library(hydra)


## Function for log-likelihood related to Jth observation
loglike <- function(clusterassign,param,data,J,n) #here J means Jth observation
{
  ################################################################
  
  ## Input: clusterassign = clustering configuration, a n by 1 vector ##
  ##        param = probability matrix, a k by k matrix ##
  ##        data = the adjacency matrix, a n by n matrix ##
  ##        J = observation index ##
  ##        n = number of observations ##
  
  ## Output: log-likelihood related to Jth observation ##
  
  #################################################################
  clustersize = max(clusterassign)
  param = as.matrix(param)
  
  if (J==1) {result2 = 0
  for (ii in c((J+1):n))
  {
    result2 = result2 + data[J,ii]*log(param[clusterassign[J],clusterassign[ii]])+(1-data[J,ii])*log(1-param[clusterassign[J],clusterassign[ii]])
  }
  output = sum(result2)} else if (J==n){
    result = 0
    for (ii in c(1:(J-1)))
    {
      result = result + data[ii,J]*log(param[clusterassign[ii],clusterassign[J]])+(1-data[ii,J])*log(1-param[clusterassign[ii],clusterassign[J]])
    }
    output = sum(result)
  } else {
    result = 0
    for (ii in c(1:(J-1)))
    {
      result = result + data[ii,J]*log(param[clusterassign[ii],clusterassign[J]])+(1-data[ii,J])*log(1-param[clusterassign[ii],clusterassign[J]])
    }
    
    result2 = 0
    for (ii in c((J+1):n))
      
    {
      result2 = result2 + data[J,ii]*log(param[clusterassign[J],clusterassign[ii]])+(1-data[J,ii])*log(1-param[clusterassign[J],clusterassign[ii]])
    }
    output = sum(result)+sum(result2)}
  output
}

#function for getting m(Aj)
logmargs <- function(clusterassign,data,J,beta.a,beta.b) #here J means Jth observation
{
  ################################################################
  
  ## Input: clusterassign = clustering configuration, a n by 1 vector ##
  ##        data = the adjacency matrix, a n by n matrix ##
  ##        J = observation index ##
  ##        n = number of observations ##
  ##        beta.a, beta.b = hyperparameters for the prior on elements in Q matrix in Beta distribution ##
  
  ## Output: m(A_j) in Algotithm 1 (collapsed sampler for MFM-SBM) ##
  
  #################################################################
  clustersize = max(clusterassign)-1
  result = NULL
  for (ii in 1:clustersize)
  {
    sumA =  sum(data[J,which(clusterassign==ii)[which(clusterassign==ii)>J]]) + sum(data[which(clusterassign==ii)[which(clusterassign==ii)<J],J])
    S = length(which(clusterassign==ii)[which(clusterassign==ii)>J]) + length(which(clusterassign==ii)[which(clusterassign==ii)<J])
    result[ii] = lbeta(sumA+beta.a,S-sumA+beta.b)-lbeta(beta.a,beta.b)
  }
  sum(result)
}

## function for Collapsed sampler for MFM-SBM (main algorithm)
CDMFM_new <- function(data, data1, niterations, beta.a, beta.b, GAMMA, LAMBDA, initNClusters)
{
  ## Model: A_{ij}|z,Q \sim Ber(Q_{z_i,z_j}) ##
  ##        Q_{rs} \sim Beta(beta.a,beta.b), r,s = 1,...,k ##
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(GAMMA,...,GAMMA) ##
  ##        k-1 \sim possion(1) ##
  
  
  
  ################################################################
  
  ## Input: data = the adjacency matrix, a n by n matrix ##
  ##        data1 = the upper traiangle for the adjacency matrix, a n by n matrix ##
  ##        niterations = the total number of iterations in MFM-SBM ##
  ##        beta.a, beta.b = hyperparameters for the prior on elements in Q matrix in Beta distribution ##
  ##        GAMMA = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        LAMBDA = the parameter for Poisson distrition ##
  ##        initNClusters = the initial number of clusters ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n by 1 vector##
  ##         Qout = probability matrix, a k by k matrix ##
  
  #################################################################
  n = dim(data)[1]
  
  #precomputation for prespecified coefficient VN -> THIS IS THE EXP OF THAT QUANTITY
  lambda <- LAMBDA
  gamma <- GAMMA
  N=n ## n is the number of oberservations
  VN<-0
  tmax = n+10
  for (t in 1:tmax)
  {
    r = log(0)
    for (k in t:500)
    {
      b = sum(log((k-t+1):k))-sum(log((k*gamma):(k*gamma+N-1))) + dpois(k-1, lambda, log = TRUE)
      m = max(b,r)
      r = log(exp(r-m) + exp(b-m)) + m
    }
    VN[t] = r
  }
  # initialization of clustering configuration
  clusterAssign <- c(sample(1:initNClusters, size = initNClusters, replace = FALSE),
                     sample(1:initNClusters, size = n-initNClusters, replace = TRUE))
  
  Q<-matrix(0, initNClusters,initNClusters)
  for (i in 1:initNClusters){
    for (j in i:initNClusters){
      Q[i,j] = rbeta(1,beta.a,beta.b)
      Q[j,i] = Q[i,j]
    }
  }
  History <- vector("list", niterations)
  
  ##start Gibb's sampling
  for (iter in 1:niterations)
  {
    ## update z ##
    clusterSizes = table(as.factor(clusterAssign))
    nClusters = length(clusterSizes)
    for (i in 1:n)
    { #determine whether ith component is a singleton 
      cur.cluster.i = clusterAssign[i]
      if (clusterSizes[clusterAssign[i]] > 1){
        # not a singleton, have |C|+1 choices
        c.counts.noi = clusterSizes  #c.counts.noi corresponds to |C|
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          (GAMMA+c.counts.noi[x])*exp(loglike(clusterAssign_temp,Q,data,i,n))
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-GAMMA*exp(logmargs(clusterAssign_1,data,i,beta.a,beta.b))*exp(VN[nClusters+1]-VN[nClusters])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        clusterAssign[i] <- cluster.i
        
        if (cluster.i > nClusters)
        {
          QQ = matrix(0,nClusters+1,nClusters+1)
          QQ[1:nClusters,1:nClusters] = Q
          QQ[nClusters+1,1:(nClusters+1)] = rbeta(nClusters+1,beta.a,beta.b)
          QQ[1:(nClusters+1),nClusters+1] = QQ[nClusters+1,1:(nClusters+1)]
          Q = QQ
          clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
          nClusters <- length(clusterSizes)} else
          {Q = Q
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes)}
      } else {
        # a singleton, have |C| choices
        c.counts.noi = clusterSizes
        c.counts.noi[clusterAssign[i]] = c.counts.noi[clusterAssign[i]] - 1 - GAMMA# can offset the gamma adding later
        #finding the probs for sampling process
        clusterProbs = sapply(1:nClusters, function(x) {
          clusterAssign_temp = clusterAssign
          clusterAssign_temp[i] = x
          (GAMMA+c.counts.noi[x])*exp(loglike(clusterAssign_temp,Q,data,i,n))
        })
        clusterAssign_1 = clusterAssign
        clusterAssign_1[i] = nClusters+1
        clusterProbs[nClusters+1]<-GAMMA*exp(logmargs(clusterAssign_1,data,i,beta.a,beta.b))*exp(VN[nClusters]-VN[nClusters-1])
        #choose the cluster number for ith observation
        cluster.i <- sample(1:(nClusters+1), size = 1,
                            prob = clusterProbs)
        # remove the empty cluster
        if (cluster.i > nClusters)
        {      clusterAssign[i] <- cur.cluster.i #put the new cluster in the place of the only singleten one
        clusterSizes <- table(as.factor(clusterAssign)) # sorts according to labels
        } else
        {      
          clusterAssign[i] <- cluster.i
          clusterAssign <- ifelse(clusterAssign > cur.cluster.i, clusterAssign-1, clusterAssign) # to delete the previous group index
          clusterSizes <- table(as.factor(clusterAssign))
          nClusters <- length(clusterSizes) 
          if (nClusters > 1) {Q = Q[-cur.cluster.i,][,-cur.cluster.i]} else {Q = Q[-cur.cluster.i,][-cur.cluster.i]}}
      }
    }
    # end for loop over subjects i
    ## update Q ##
    Q = matrix(0, nClusters,nClusters)
    AA = matrix(0,nClusters,nClusters)
    NN = matrix(0,nClusters,nClusters)
    for (r in 1:nClusters){
      for (s in r:nClusters)
      {
        AA[r,s] = sum(data1[clusterAssign==r,clusterAssign==s]) + sum(data1[clusterAssign==s,clusterAssign==r]) - (r==s)*sum(data1[clusterAssign==s,clusterAssign==r])
        med = matrix(0,n,n)
        med[which(clusterAssign==r),which(clusterAssign==s)] = 1
        med1 = matrix(0,n,n)
        med1[which(clusterAssign==s),which(clusterAssign==r)] = 1
        NN[r,s] = sum(med*lower.tri(med)) + sum(med1*lower.tri(med1))-(r==s)*sum(med1*lower.tri(med1))
        Q[r,s] = rbeta(1,AA[r,s]+beta.a,NN[r,s]-AA[r,s]+beta.b)
        Q[s,r] = Q[r,s]
      }
    }
    History[[iter]] <- list(zout = clusterAssign,Qout = Q)
    print(table(clusterAssign))
    #cat(" iteration:", iter,"\n",clusterAssign,"\n")
  }# for loop over iterations
  
  list(Iterates = History)
}

## Dahl's method to summarize the samples from the MCMC
getDahl <- function(MFMfit, burn)
{
  ################################################################
  
  ## Input: MFMfit = the result from CDMFM_new ## 
  ##        burn = the number of burn-in iterations ##
  
  ## Output: 
  ##         zout = estimated clustering configuration, a n by 1 vector##
  ##         Qout = estimated probability matrix, a k by k matrix ##
  
  #################################################################
  iters <- MFMfit$Iterates[-(1:burn)]
  n <- length(iters[[1]][[1]])
  niters <- length(iters)
  membershipMatrices <- lapply(iters, function(x){
    clusterAssign <- x[[1]]
    outer(clusterAssign, clusterAssign, FUN = "==")
  })
  membershipAverage <- Reduce("+", membershipMatrices)/niters
  SqError <- sapply(membershipMatrices, function(x, av) sum((x - av)^2),
                    av = membershipAverage)
  DahlIndex <- which.min(SqError)
  DahlAns <- iters[[DahlIndex]]
  attr(DahlAns, "iterIndex") <- burn + DahlIndex
  attr(DahlAns, "burnin") <- burn
  DahlAns
}



## dolphin data example
load("/home/michele/Desktop/Jasa-Article/dolphindata.RData")

## loading the data
diag(AAA) = 0 ## make it without-selfloop network
diag(A) = 0 ## make it without-selfloop networkk

fit_dol = CDMFM_new(data = A, data1 = AAA, niterations = 10000, beta.a = 1, beta.b = 1, GAMMA=1,LAMBDA = 1,initNClusters=9)
result_dol = getDahl(fit_dol, burn = 4000)

set.seed(3)
fit_dol = CDMFM_new(data = A, data1 = AAA, niterations = 10000, beta.a = 2, beta.b = 2, GAMMA=1,LAMBDA = 1,initNClusters=ceiling(runif(1,1,10)))
    result_dol = getDahl(fit_dol, burn = 4000)


# My function
a = b = 2
gamma_ = 1
lambda = 1
K_start = 4
set.seed(5)
logVn.Miller <- log_Vn.M(gamma_ = 1, n, n+10,lambda)
z_start       = c(sample(1:K_start, size = K_start, replace = FALSE),
                  sample(1:K_start, size = n-K_start, replace = TRUE))
fit.po = fcogibbs_po(40000, 9,A, n=ncol(A), a, b, gamma_ = 1, logVn.Miller, z0 = NULL, z_start = z_start)
hist(fit.po$num_k, prob = T)

fit.gn = fcogibbs_gn(40000, 9,A, n=ncol(A), a, b, gamma_ = 1, gamma.gn = 0.1, z0 = NULL, z_start = z_start)
hist(fit.gn$num_k[-(1:10000)], prob = T)


## estimated configuration using Dahl's method, choosing first 100 iterations in MCMC as burn-in
result_dol = getDahl(fit.po, burn = 4000)

