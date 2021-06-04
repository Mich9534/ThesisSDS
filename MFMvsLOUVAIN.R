rm(list = ls())
setwd('/home/michele/Desktop/Small_code_R')
#---------------------------------------------
# MFM-SBM ALGORITHM VS LOUVAIN ALGORITHM
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
#---------------------------------------------
n      = 100
K      = 2                 # we can change as in Geng et al
p      = 0.24    
q      = 0.1     
M      = 600     
BI     = 300   
N      = 100
K.l    = rep(NA, N)        # Estimation of K through louvain algorithm              
K.g    = rep(NA, N)        # Estimation of K through the gibbs sampler
VI.l   = rep(NA, N)        # Variation Information for Louvain
VI.g   = rep(NA, N)        # Variation Information for Gibbs
RI.l   = rep(NA, N)        # Rand Index for Louvain
RI.g   = rep(NA, N)        # Rand Index for Gibbs
lambda = 1
gamma_ = a = b = 1
K_start = 9

find_mode = function(vect){
  WHO.max <- as.numeric(names(which(table(vect) == max(table(vect)))))
  n.max   <- length(WHO.max)
  if(n.max != 1){
    u = sample(1:n.max,1)
    MODE = WHO.max[u]
  } else{MODE = WHO.max}
  return(MODE)
}


for(i in 1:N){
  print(paste("iterazione:", i))
  # Gibss sampler
  print(paste("Turno di MFM-SBM"))
  data_gen <- fgenera_data(n, K, p, q, type_network = 'balanced')
  A        <- data_gen$A0
  Z0       <- data_gen$z0
  z_start  <- c(sample(1:K_start, size = K_start, replace = FALSE),
               sample(1:K_start, size = n-K_start, replace = TRUE))
  logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
  fit.po        = fcogibbs_po(M, K_start, A, n, a, b, gamma_, logVn.Miller, z_start)
  K.g[i]           = find_mode(fit.po$num_k[-c(1:BI)])
  # K.g = find_mode(fit.po$num_k)
  c_Z_PO <- pr_cc(fit.po$z_post[,(BI+1):M])
  memb_Z_PO_VI <- minVI(c_Z_PO,method="avg",max.k=20) # from package mcclust.ext
  memb_Z_PO <- memb_Z_PO_VI$cl
  VI.g[i] <- compare(memb_Z_PO, Z0, method = "vi")
  RI.g[i] <- compare(memb_Z_PO, Z0, method = "rand")
  
  
  # Louvain algorithm
  print(paste("Turno di Louvain"))
  net    <- graph.adjacency(A, mode=c("undirected"), weighted=NULL, diag=FALSE)   # transform the adjacency matrix into an igraph object
  Louv   <- cluster_louvain(net)$membership                                       # point estimate
  K.l[i] <- length(table(Louv))                                                   # estimated K
  VI.l[i] <- compare(Louv, Z0, method = "vi")
  RI.l[i] <- compare(Louv, Z0, method = "rand")
  
}

df = data.frame(K.l, K.g)
df2 = data.frame(x = rep(2, N), y = rep(2,N))
hist(K.l)


p1 = ggplot(df, aes(x=K.l)) + geom_bar()+theme_bw()+ylim(0,100)+xlim(2,7.5)+ggtitle("LOUVAIN") +
  xlab("estimated # of communities")+theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p2 = ggplot(df, aes(x=K.g)) + geom_bar()+theme_bw()+ylim(0,100)+xlim(1,6.5)+ggtitle("MFM-SBM") +
  xlab("estimated # of communities")+theme(plot.title = element_text(hjust = 0.5, face = "bold"))

p3 = ggplot(df2, aes(x=x)) + geom_bar()+theme_bw()+ylim(0,100)+xlim(1,5)+ggtitle("LOUVAIN") +
  xlab("estimated # of communities")+theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p4 = ggplot(df2, aes(x=y)) + geom_bar()+theme_bw()+ylim(0,100)+xlim(1,5)+ggtitle("MFM-SBM") +
  xlab("estimated # of communities")+theme(plot.title = element_text(hjust = 0.5, face = "bold"))


png("MvL_K2_5.png", width = 718, height = 513)
grid.arrange(p3, p4, ncol=2)
dev.off()

mean(RI.g)
table(K.l)/100
