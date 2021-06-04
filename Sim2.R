# SENSITIVITY ANALYSIS - SIM2 DATA.R
rm(list = ls())
load("~/Desktop/Small_code_R/Sim2.RData")
# SETTING:
# n = 100
# K = 2
# N = 100
# M = 600
# BI = 300
# UNBALANCED NETWORK
# K_start = 9
# p = 0.24
# q = 0.1
# a = b = gamma = 1
# 0.1 0.3 0.5 0.7 0.9 -> gnedin  par
#  1   3   5   7   9  -> poisson par

library(extraDistr)
library(ggplot2)
theme_set(theme_bw())
N = 100
#-------------------------------------------------------
# K ESTIMATE COMPARISON GNEDIN - TR POISSON
#-------------------------------------------------------
c1_k.gn = cbind(s.a.$k_es.gn[,1], rep(0.1, N))
c2_k.gn = cbind(s.a.$k_es.gn[,2], rep(0.3, N))
c3_k.gn = cbind(s.a.$k_es.gn[,3], rep(0.5, N))
c4_k.gn = cbind(s.a.$k_es.gn[,4], rep(0.7, N))
c5_k.gn = cbind(s.a.$k_es.gn[,5], rep(0.9, N))

c.whole_k.gn = rbind(c1_k.gn,c2_k.gn,c3_k.gn,c4_k.gn,c5_k.gn)
df.k_es.gn = data.frame(c.whole_k.gn)
df.k_es.gn$X2 = as.factor(df.k_es.gn$X2)

p_k.gn <- ggplot(df.k_es.gn, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("Gnedin's gamma")+ylab(expression(paste(hat(K), " (K = 2)")))
#p_k.gn

c1_k.po = cbind(s.a.$k_es.po[,1], rep(1, N))
c2_k.po = cbind(s.a.$k_es.po[,2], rep(3, N))
c3_k.po = cbind(s.a.$k_es.po[,3], rep(5, N))
c4_k.po = cbind(s.a.$k_es.po[,4], rep(7, N))
c5_k.po = cbind(s.a.$k_es.po[,5], rep(9, N))
c.whole_k.po = rbind(c1_k.po,c2_k.po,c3_k.po,c4_k.po,c5_k.po)
df.k_es.po = data.frame(c.whole_k.po)
df.k_es.po$X2 = as.factor(df.k_es.po$X2)

p_k.po <- ggplot(df.k_es.po, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("tr Poisson's lambda")+ylab(expression(paste(hat(K), " (K = 2)")))
#p_k.po

require(gridExtra)
grid.arrange(p_k.gn, p_k.po, nrow=2)


#-------------------------------------------------------
# Post Prob ESTIMATE COMPARISON GNEDIN - TR POISSON
#-------------------------------------------------------
c1_prob.gn = cbind(s.a.$prob.gn[,1], rep(0.1, N))
c2_prob.gn = cbind(s.a.$prob.gn[,2], rep(0.3, N))
c3_prob.gn = cbind(s.a.$prob.gn[,3], rep(0.5, N))
c4_prob.gn = cbind(s.a.$prob.gn[,4], rep(0.7, N))
c5_prob.gn = cbind(s.a.$prob.gn[,5], rep(0.9, N))
c.whole_prob.gn = rbind(c1_prob.gn,c2_prob.gn,c3_prob.gn,c4_prob.gn,c5_prob.gn)
df.prob.gn = data.frame(c.whole_prob.gn)
df.prob.gn$X2 = as.factor(df.prob.gn$X2)

p_prob.gn <- ggplot(df.prob.gn, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("Gnedin's gamma")+ylab(expression(paste("P(",hat(K), "= 2 | A)")))
#p_prob.gn

c1_prob.po = cbind(s.a.$prob.po[,1], rep(1, N))
c2_prob.po = cbind(s.a.$prob.po[,2], rep(3, N))
c3_prob.po = cbind(s.a.$prob.po[,3], rep(5, N))
c4_prob.po = cbind(s.a.$prob.po[,4], rep(7, N))
c5_prob.po = cbind(s.a.$prob.po[,5], rep(9, N))
c.whole_prob.po = rbind(c1_prob.po,c2_prob.po,c3_prob.po,c4_prob.po,c5_prob.po)
df.prob.po = data.frame(c.whole_prob.po)
df.prob.po$X2 = as.factor(df.prob.po$X2)

p_prob.po <- ggplot(df.prob.po, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("tr Poisson's lambda")+ylab(expression(paste("P(",hat(K), "= 2 | A)")))
#p_prob.po

require(gridExtra)
grid.arrange(p_prob.gn, p_prob.po, nrow=2)

#-------------------------------------------------------
# EXECUTION TIME ESTIMATE COMPARISON GNEDIN - TR POISSON
#-------------------------------------------------------
c1_TIME.gn = cbind(s.a.$TIME.gn[,1], rep(0.1, N))
c2_TIME.gn = cbind(s.a.$TIME.gn[,2], rep(0.3, N))
c3_TIME.gn = cbind(s.a.$TIME.gn[,3], rep(0.5, N))
c4_TIME.gn = cbind(s.a.$TIME.gn[,4], rep(0.7, N))
c5_TIME.gn = cbind(s.a.$TIME.gn[,5], rep(0.9, N))
c.whole_TIME.gn = rbind(c1_TIME.gn,c2_TIME.gn,c3_TIME.gn,c4_TIME.gn,c5_TIME.gn)
df.TIME.gn = data.frame(c.whole_TIME.gn)
df.TIME.gn$X2 = as.factor(df.TIME.gn$X2)

p_TIME.gn <- ggplot(df.TIME.gn, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("Gnedin's gamma")+ylab("Execution Time (sec)")
#p_TIME.gn

c1_TIME.po = cbind(s.a.$TIME.po[,1], rep(1, N))
c2_TIME.po = cbind(s.a.$TIME.po[,2], rep(3, N))
c3_TIME.po = cbind(s.a.$TIME.po[,3], rep(5, N))
c4_TIME.po = cbind(s.a.$TIME.po[,4], rep(7, N))
c5_TIME.po = cbind(s.a.$TIME.po[,5], rep(9, N))
c.whole_TIME.po = rbind(c1_TIME.po,c2_TIME.po,c3_TIME.po,c4_TIME.po,c5_TIME.po)
df.TIME.po = data.frame(c.whole_TIME.po)
df.TIME.po$X2 = as.factor(df.TIME.po$X2)

p_TIME.po <- ggplot(df.TIME.po, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("tr Poisson's lambda")+ylab("Execution Time (sec)")
#p_TIME.po

require(gridExtra)
grid.arrange(p_TIME.gn, p_TIME.po, nrow=2)

#-------------------------------------------------------
# RAND INDEX ESTIMATE COMPARISON GNEDIN - TR POISSON
#-------------------------------------------------------
c1_RaIn.gn = cbind(s.a.$RaIn.gn[,1], rep(0.1, N))
c2_RaIn.gn = cbind(s.a.$RaIn.gn[,2], rep(0.3, N))
c3_RaIn.gn = cbind(s.a.$RaIn.gn[,3], rep(0.5, N))
c4_RaIn.gn = cbind(s.a.$RaIn.gn[,4], rep(0.7, N))
c5_RaIn.gn = cbind(s.a.$RaIn.gn[,5], rep(0.9, N))
c.whole_RaIn.gn = rbind(c1_RaIn.gn,c2_RaIn.gn,c3_RaIn.gn,c4_RaIn.gn,c5_RaIn.gn)
df.RaIn.gn = data.frame(c.whole_RaIn.gn)
df.RaIn.gn$X2 = as.factor(df.RaIn.gn$X2)

p_RaIn.gn <- ggplot(df.RaIn.gn, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("Gnedin's gamma")+ylab("RI Index on point est.")+ylim(0.5,1)

#p_RaIn.gn

c1_RaIn.po = cbind(s.a.$RaIn.po[,1], rep(1, N))
c2_RaIn.po = cbind(s.a.$RaIn.po[,2], rep(3, N))
c3_RaIn.po = cbind(s.a.$RaIn.po[,3], rep(5, N))
c4_RaIn.po = cbind(s.a.$RaIn.po[,4], rep(7, N))
c5_RaIn.po = cbind(s.a.$RaIn.po[,5], rep(9, N))
c.whole_RaIn.po = rbind(c1_RaIn.po,c2_RaIn.po,c3_RaIn.po,c4_RaIn.po,c5_RaIn.po)
df.RaIn.po = data.frame(c.whole_RaIn.po)
df.RaIn.po$X2 = as.factor(df.RaIn.po$X2)

p_RaIn.po <- ggplot(df.RaIn.po, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("tr Poisson's lambda")+ylab("RI Index on point est.")+ylim(0.5,1)
#p_RaIn.po

require(gridExtra)
grid.arrange(p_RaIn.gn, p_RaIn.po, nrow=2)

#-------------------------------------------------------
# VI INDEX ESTIMATE COMPARISON GNEDIN - TR POISSON
#-------------------------------------------------------
c1_VaIn.gn = cbind(s.a.$peVI.gn[,1], rep(0.1, N))
c2_VaIn.gn = cbind(s.a.$peVI.gn[,2], rep(0.3, N))
c3_VaIn.gn = cbind(s.a.$peVI.gn[,3], rep(0.5, N))
c4_VaIn.gn = cbind(s.a.$peVI.gn[,4], rep(0.7, N))
c5_VaIn.gn = cbind(s.a.$peVI.gn[,5], rep(0.9, N))
c.whole_VaIn.gn = rbind(c1_VaIn.gn,c2_VaIn.gn,c3_VaIn.gn,c4_VaIn.gn,c5_VaIn.gn)
df.VaIn.gn = data.frame(c.whole_VaIn.gn)
df.VaIn.gn$X2 = as.factor(df.VaIn.gn$X2)

p_VaIn.gn <- ggplot(df.VaIn.gn, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("Gnedin's gamma")+ylab("VI Index on point est.")+ylim(0, 0.7)
#p_RaIn.gn

c1_VaIn.po = cbind(s.a.$peVI.po[,1], rep(1, N))
c2_VaIn.po = cbind(s.a.$peVI.po[,2], rep(3, N))
c3_VaIn.po = cbind(s.a.$peVI.po[,3], rep(5, N))
c4_VaIn.po = cbind(s.a.$peVI.po[,4], rep(7, N))
c5_VaIn.po = cbind(s.a.$peVI.po[,5], rep(9, N))
c.whole_VaIn.po = rbind(c1_VaIn.po,c2_VaIn.po,c3_VaIn.po,c4_VaIn.po,c5_VaIn.po)
df.VaIn.po = data.frame(c.whole_VaIn.po)
df.VaIn.po$X2 = as.factor(df.VaIn.po$X2)

p_VaIn.po <- ggplot(df.VaIn.po, aes(x=X2, y=X1)) + 
  geom_boxplot()+xlab("tr Poisson's lambda")+ylab("VI Index on point est.")+ylim(0, 0.7)
#p_RaIn.po

require(gridExtra)
grid.arrange(p_VaIn.gn, p_VaIn.po, nrow=2)

