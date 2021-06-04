rm(list = ls())
setwd('/home/michele/Desktop/Small_code_R')
#---------------------------------------------
# SENSITIVITY ANALYSIS OF THE TWO ALGORITHM
#---------------------------------------------
source("cogibbs_gn.R")
source("cogibbs_po.R")
source("fgenera_data.R")
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

n = 100
K = 2
p = 0.5
q = 0.1
K_start = 9
gamma_ = 1
lambda = 1
a = b = 1
M = 600
N = 100
BI = 300

find_mode = function(vect){
  WHO.max <- as.numeric(names(which(table(vect) == max(table(vect)))))
  n.max   <- length(WHO.max)
  if(n.max != 1){
    u = sample(1:n.max,1)
    MODE = WHO.max[u]
  } else{MODE = WHO.max}
  return(MODE)
}

n.vect = seq(20, 100, 20)
prob.po = matrix(NA, nrow = length(n.vect), ncol = N)
for (i in 1:length(n.vect)) {
  print(paste("n:", n.vect[i]))
  n = n.vect[i]
  for(ii in 1:N){
    print(paste("ii:", ii))
    data_gen <- fgenera_data(n, K, p, q, type_network = 'unbalanced')
    A = data_gen$A0
    Z0  <- data_gen$z0
    z_start <- c(sample(1:K_start, size = K_start, replace = FALSE),
                 sample(1:K_start, size = n-K_start, replace = TRUE))
    
    logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
    fit.po        = fcogibbs_po(M, K_start, A, n, a, b, gamma_, logVn.Miller, z_start)
    #print(fit.po$num_k)
    prob.po[i,ii]  = mean(fit.po$num_k[-c(1:BI)] == K)
  }
}

p.avg = apply(prob.po, 1, mean)
plot(n.vect, p.avg, type = "l")
# data_gen <- fgenera_data(n, K, p, q, type_network = 'unbalanced')
# A = data_gen$A0
# Z0  <- data_gen$z0
# z_start <- c(sample(1:K_start, size = K_start, replace = FALSE),
#              sample(1:K_start, size = n-K_start, replace = TRUE))
# 
# logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
# fit.po        = fcogibbs_po(M, K_start, A, n, a, b, gamma_, logVn.Miller, z_start)
# prob.po[i,ii]  = mean(fit.po$num_k[-c(1:BI)] == K)
