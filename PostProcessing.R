# Post-processing graph/probability p(k|t) (3.7) formula
rm(list = ls())
setwd('/home/michele/Desktop/Small_code_R')
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
library('extraDistr')
source('Vn_Miller.R')
source('fgenera_data.R')
source('cogibbs_po.R')
source('vec2mat.R')
source('pr_cc.R')
source('edge_est.R')
source('lVnGn.R')

# Compute p[k,t]=p(k|t) under the given MFM parameters, for k=1:upto_k, t=1:upto_t.
post_k <- function(gamma_,n,upto_k,upto_t){
  if(upto_t>n){print("p(k|t) is undefined for t>n.")}else{
    log_v <- log_Vn.M(gamma_, n, upto_t, 1)
    p = matrix(rep(0, upto_k*upto_t), upto_k, upto_t)
    for(k in 1:upto_k){
      for (t in 1:min(k, upto_t)) {
        b = lgamma(k+1)-lgamma(k-t+1)-lgamma(k*gamma_+n)+lgamma(k*gamma_)+dtpois(k, lambda = 1, a = 0,b = Inf, log = T)
        p[k, t] = exp(b-log_v[t])
        
      }
    }
  return(p)}
}


find_mode = function(vect){
  WHO.max <- as.numeric(names(which(table(vect) == max(table(vect)))))
  n.max   <- length(WHO.max)
  if(n.max != 1){
    u = sample(1:n.max,1)
    MODE = WHO.max[u]
  } else{MODE = WHO.max}
  return(MODE)
}

#-----------------------------------------------------------------------------------
# # Compute posterior on k (# of components)
# function k_posterior(result; upto=result.options.t_max)
#   o = result.options
#   @assert(o.model_type=="MFM", "The posterior on k is not defined for the DPM.")
#   lpk = eval(Meta.parse(o.log_pk))
#   log_pk_fn(k) = Base.invokelatest(lpk,k)
#   p_kt = MFM.p_kt(log_pk_fn,o.gamma,o.n,upto,o.t_max)
#   return p_kt*t_posterior(result)
# end
#         NOTE: post processing should be made on the single chain!!!
#-----------------------------------------------------------------------------------
gamma_ <- 1
lambda = 1
a = b = 1
p = 0.5
q = 0.1
K = 4
M = 300
K_start = 9
BI = 100
n = 150

K_tot <- rep(NA, 5)
K_chains <- rep(NA, 5)
for (iter in 1:5) {
  print(paste("dataset N.:", iter))
  data_gen <- fgenera_data(n, K, p, q, type_network = "Balanced")    # Let's generate a new dataset with the same specifics
  A <- data_gen$A0
  Z0  <- data_gen$z0
  z_start <- c(sample(1:K_start, size = K_start, replace = FALSE),
               sample(1:K_start, size = n-K_start, replace = TRUE))
  
  logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
  fit.po        = fcogibbs_po(M, K_start, A, n, a, b, gamma_, logVn.Miller, z_start)
  # fit = CDMFM_new(data = A, data1 = AAA, niterations = M, beta.a = 1, beta.b = 1, GAMMA=1, LAMBDA = 1, initNClusters = K_start)
  K_chains[iter] <- find_mode(fit.po$num_k[-c(1:BI)])
  
  # c_Z_PO <- pr_cc(fit.po$z_post[,(BI+1):M])
  # memb_Z_PO_VI <- minVI(c_Z_PO,method="avg",max.k=20) # from package mcclust.ext
  # memb_Z_PO <- memb_Z_PO_VI$cl
  # table(memb_Z_PO)
  # table(fit.po$num_k)
}

max.t  <- n
p_kt   <- post_k(gamma_, n, upto_k = max.t, upto_t = max.t)
p_p.ch <- rep(NA, max.t) 
fp_p.ch<- function(p_p.ch){
  for(i in 1:max.t) {
  p_p.ch[i] = mean(K_chains == i)
  }
  return(p_p.ch)
}
p_p.ch = fp_p.ch(p_p.ch)
x = 1:length(p_p.ch)
df= data.frame(p_p.ch, x)

p.gg = ggplot(df, aes(y = p_p.ch,x =x))+geom_point(type = "l", col = "#F26D6D")+xlim(1,8)+ylim(0,1)+geom_line(alpha = 0.5, col = "#F26D6D")+
  xlab("k (number of components)")+ylab("p(k | data)")+theme_bw()

p.gg2 = p.gg+geom_point(data = df, type = "l", col = "#00BFC4")+geom_line(data = df, alpha = 0.5, col = "#00BFC4")

p.gg3 = p.gg2+geom_point(data = df, type = "l", col = "#9BDB4B")+geom_line(data = df, alpha = 0.5, col = "#9BDB4B")
p.gg3+geom_segment(aes( x = 6, y = 0.8, xend = 6.25, yend = 0.8), col = "#00BFC4")+
  annotate(geom = "text", x = 6.3, y = 0.8, label = "n = 100",size = 3,hjust = "left")+
  geom_segment(aes( x = 6, y = 0.7, xend = 6.25, yend = 0.7), col = "#9BDB4B")+
  annotate(geom = "text", x = 6.3, y = 0.7, label = "n = 150",size = 3,hjust = "left")+
  geom_segment(aes( x = 6, y = 0.9, xend = 6.25, yend = 0.9), col = "#F26D6D")+
  annotate(geom = "text", x = 6.3, y = 0.9, label = "n = 50",size = 3,hjust = "left")
  
  
# k.pos <- p_kt %*% p_p.ch
# plot(k.pos, type = "l", ylim = c(0,1),xlim = c(1, 8))
# points(k.pos)
