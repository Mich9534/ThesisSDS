# Rand Index Evolution (mean of rand index of the partition estimated by means VI and the true partition
rm(list = ls())
library('ggplot2')
load("~/Desktop/Small_code_R/CONV_RI_K2_BAL.RData")
# SETTING:
# WE USE POISSON DISTRIBUTED PRIOR, TRUNCATED AT ZERO
# WE FOR EVERY n WE GENERATE N NEW DATASET, FOR EACH OF THEM WE COMPUTE AFTER EACH MCMC
# THE POINT ESTIMATION AS WRITTEN IN LEGRAMANTI DURANTE RIGON, AND WE TOOK THE AVERAGE FOR EACH n 
# K      = 2
# a = b  = 1
# gamma_ = 1 
# lambda = 1
# p      = 0.1
# q      = 0.24            # vague structure
# n      = seq(20,250,40)
# N      = 100
# K_st   = 9
# M      = 10000
# burn_in     = 4000

plot(n,t, type = "l")
D = data.frame(cbind(t, n))
ggplot(data = D, aes(x = n, y = t))+geom_point()+geom_line()+theme_light()+ylab('Rand Index')+ylim(0,1)
