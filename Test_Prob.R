rm(list = ls())
setwd('/home/michele/Desktop/Small_code_R')
#----------------------------------
# Libraries to load
#----------------------------------
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
library(mcclust)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(coda)
library(dummies)
library(randnet) 
# Setting
#----------------------------------
K      = 2
a = b  = 1
gamma_ = 1 
lambda = 1
p      = 0.1
q      = 0.5            # vague structure
n      = seq(10,100,15)
N      = 1
K_st   = 9
M      = 10000
burn_in     = 4000

test_cv = function(K, a, b, gamma_, lambda, p, q, n, N, K_start, M, burn_in){
  pr.avg = rep(NA, length(n))
  for(num in 1:length(n)){
    # num = 1
    print(n[num])
    PR.i = rep(NA, N)
    for(i in 1:N){
      # i = 1
      print(i)
      
      DATA = fgenera_data(n[num], K, p, q, 'balanced')
      A = DATA$A0
      Z0 = DATA$z0
      z_start       = c(sample(1:K_start, size = K_start, replace = FALSE),
                        sample(1:K_start, size = n[num]-K_start, replace = TRUE))
      logVn.Miller = log_Vn.M(gamma_, n[num], n[num]+10, lambda)
      fit_po = fcogibbs_po(M, K_st, A, n[num], a, b, gamma_, logVn.Miller, z0 = NULL, z_start)
      z_post_po     = fit_po$z_post
      PR.i[i] = mean(z_post_po[-c(1:burn_in)] == K)
      
    }
    pr.avg[num] = mean(PR.i)
  }
  return(pr.avg)
}

pr_2 = test_cv(K, a, b, gamma_, lambda, p, q, n, N, K_st, M, burn_in)
D_2 = data.frame(cbind(pr_2, n))
ggplot(data = D_2, aes(x = n, y = pr_2))+geom_point()+geom_line()+theme_light()+ylab('Rand Index')
