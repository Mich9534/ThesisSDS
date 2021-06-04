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

sens_anal <- function(n, N, M, K, K_start,  p, q, a, b,gamma_, BI){
  par_gn           <-  seq(0.1, 0.9, 0.2)
  par_po           <-  seq(1,10, 2)
  n.par = length(par_po)
  # res_gn           <-  data.frame(matrix(rep(NA,  10*4), ncol = 4, nrow = 10), row.names = par_gn)
  # colnames(res_gn) <- c("pr", "K.hat", "Time", "RI")
  # res_po           <-  data.frame(matrix(rep(NA,  10*4), ncol = 4, nrow = 10), row.names = par_po)
  # colnames(res_po) <- c("pr", "K.hat", "Time", "RI")
  
  TIME.po  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N)
  prob.po  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N)
  k_es.po  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N)
  RaIn.po  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N)
  peVI.po  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N) # vector of VI point est. and Z0 for each MCMC
  
  TIME.gn  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N)
  prob.gn  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N)
  k_es.gn  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N)
  RaIn.gn  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N)
  peVI.gn  <- matrix(rep(NA, N*n.par), ncol = n.par, nrow = N) # vector of VI point est. and Z0 for each MCMC
  
  find_mode = function(vect){
    WHO.max <- as.numeric(names(which(table(vect) == max(table(vect)))))
    n.max   <- length(WHO.max)
    if(n.max != 1){
      u = sample(1:n.max,1)
      MODE = WHO.max[u]
    } else{MODE = WHO.max}
    return(MODE)
  }
  
  mean.na   = function(vect.na){
    vect = na.omit(vect.na)
    na.count = length(vect.na) - length(vect)
    if(na.count/length(vect.na) < 0.5){return(mean(vect))}else{return(NA)}
  }
  
  for(i in 1:N){
    print(paste("iteration:",i))
    data_gen <- fgenera_data(n, K, p, q, type_network = 'unbalanced')
    A        <- data_gen$A0
    Z0       <- data_gen$z0
    z_start <- c(sample(1:K_start, size = K_start, replace = FALSE),
                 sample(1:K_start, size = n-K_start, replace = TRUE))
    print("ho generato i dati")
    
    for(ii in 1:length(par_po)){
      print(paste("par poisson:", par_po[ii]))
      time.start    = Sys.time()    
      logVn.Miller  = log_Vn.M(gamma_, n, n+10, par_po[ii])
      print("entro in fit.po!")
      fit.po        = fcogibbs_po(M, K_start, A, n, a, b, gamma_, logVn.Miller, z_start)
      time.end      = Sys.time()
      
      TIME.po[i,ii]  = time.end-time.start
      prob.po[i,ii]  = mean(fit.po$num_k[-c(1:BI)] == K)
      k_es.po[i,ii]  = find_mode(fit.po$num_k[-c(1:BI)])
      #if(k_es.po[i, ii] == K){RaIn.po[i, ii] = rand.index(Z0, fit.po$z_update)}
      
      c_Z_PO <- pr_cc(fit.po$z_post[,(BI+1):M])
      memb_Z_PO_VI <- minVI(c_Z_PO,method="avg",max.k=20) # from package mcclust.ext
      memb_Z_PO <- memb_Z_PO_VI$cl
      peVI.po[i,ii] <- compare(memb_Z_PO, Z0, method = "vi")
      RaIn.po[i,ii] <- compare(memb_Z_PO, Z0, method = "rand")
    }
    
    for(jj in 1:length(par_gn)){
      print(paste("par gnedin:", par_gn[jj]))
      time.start     = Sys.time()    
      print("entro in fit.gn!")
      fit.gn         = fcogibbs_gn(M, K_start, A, n, a, b, gamma_, par_gn[jj],  z_start)
      time.end       = Sys.time()
      
      TIME.gn[i,jj]  = time.end-time.start
      prob.gn[i,jj]  = mean(fit.gn$num_k[-c(1:BI)] == K)
      k_es.gn[i,jj]  = find_mode(fit.gn$num_k[-c(1:BI)])
      #if(k_es.gn[i, jj] == K){RaIn.gn[i, jj] = rand.index(Z0, fit.gn$z_update)}
      
      c_Z_GN <- pr_cc(fit.gn$z_post[,(BI+1):M])
      memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20) # from package mcclust.ext
      memb_Z_GN <- memb_Z_GN_VI$cl
      peVI.gn[i,jj] <- compare(memb_Z_GN, Z0, method = "vi")
      RaIn.gn[i,jj] <- compare(memb_Z_GN, Z0, method = "rand")
    }
  }
  # res_po[,1] = prob.po/N
  # res_po[,2] = apply(k_es.po, 2, find_mode)
  # res_po[,3] = TIME.po/N
  # res_po[,4] = apply(RaIn.po, 2, mean.na)
  # 
  # res_gn[,1] = prob.gn/N
  # res_gn[,2] = apply(k_es.gn, 2, find_mode)
  # res_gn[,3] = TIME.gn/N
  # res_gn[,4] = apply(RaIn.gn, 2, mean.na)
  results = list(TIME.gn, TIME.po, prob.gn, prob.po, k_es.gn, k_es.po, RaIn.gn, RaIn.po, peVI.gn, peVI.po)
  names(results) = c('TIME.gn', 'TIME.po', 'prob.gn', 'prob.po', 'k_es.gn', 'k_es.po', 'RaIn.gn', 'RaIn.po', 'peVI.gn', 'peVI.po')
  return(results)
}
#s.a. = function(n, N, M, K, K_start,  p, q, a, b,gamma_, BI)
s.a. = sens_anal(200, 100, 600, K = 5, K_start= 9, p = 0.24, q = 0.1, 1, 1, 1, 300)
s.a.$k_es.po
s.a.$k_es.gn
#--------------------------------------------------
# PLOTTING THE RESULTS
#--------------------------------------------------

DF1 <- s.a.$resPO
DF2 <- s.a.$resGN

par_gn  <-  seq(0.1, 0.9, 0.085)
par_po  <-  seq(1,20, 2)


ggplot()+
  geom_line(data = DF1,
                 aes(x=par_po, y=pr), color = "#598568")+
  ylim(0,1)+
  scale_x_continuous(breaks = par_po)

ggplot()+
  geom_line(data = DF2,
            aes(x=par_gn, y=pr), color = "#C74E4A")+
  ylim(0,1)+
  scale_x_continuous(breaks = par_gn)


ggplot()+
  geom_linerange(data = DF1,
            aes(x=par_po, ymax=K.hat, ymin = 0), color = "#598568")+
  scale_x_continuous(breaks = par_po)


ggplot()+
  geom_linerange(data = DF2,
                 aes(x=par_gn, ymax=K.hat, ymin = 0), color = "#C74E4A")+
  ylim(0, 3)+
  scale_x_continuous(breaks = par_gn)


ggplot()+
  geom_line(data = DF1,
                 aes(x=par_po, y=Time), color = "#598568")+
  ylim(0, 12)+
  scale_x_continuous(breaks = par_po)


ggplot()+
  geom_line(data = DF2,
            aes(x=par_gn, y=Time), color = "#C74E4A")+
  ylim(0, 12)+
  scale_x_continuous(breaks = par_gn)



ggplot()+
  geom_linerange(data = DF1,
                 aes(x=par_po, ymax=RI, ymin = 0), color = "#598568")+
  ylim(0,1)+
  scale_x_continuous(breaks = par_po)


ggplot()+
  geom_linerange(data = DF2,
                 aes(x=par_gn, ymax=RI, ymin = 0), color = "#C74E4A")+
  ylim(0, 1)+
  scale_x_continuous(breaks = par_gn)


  geom_linerange(data = DF2,
                 aes(x=x2, ymax = y2, ymin = 0), color = "#C74E4A")+
  geom_hline(yintercept = mean(s.a.$pr.po-0.01), linetype = "dashed", color = '#598568', size = 1.05)+
  geom_hline(yintercept = mean(s.a.$pr.gn+0.005), linetype = "dashed", color = "#C74E4A", size = 1.05)+
  theme(panel.background = element_rect(fill= "#D1CFD0"), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks = round(c(0,mean(s.a.$pr.po)-0.0065, mean(s.a.$pr.gn)+0.009, 1),2))+
  xlab('iteration')+
  ylab(expression(paste("Mean of  ", Pi("k = K | A"))))
  
#...Da qua....
DF1 <- data.frame(x1=1:100, y1=s.a.$K.po)
DF2 <- data.frame(x2=(1:100)+0.3, y2=s.a.$K.gn)


ggplot()+
  geom_linerange(data = DF1,
                 aes(x=x1, ymax=y1, ymin = 0 ), color = "#598568")+
  geom_linerange(data = DF2,
                 aes(x=x2, ymax = y2, ymin = 0), color = "#C74E4A")+
  geom_hline(yintercept = mean(s.a.$pr.po-0.01), linetype = "dashed", color = '#598568', size = 1.05)+
  geom_hline(yintercept = mean(s.a.$pr.gn+0.005), linetype = "dashed", color = "#C74E4A", size = 1.05)+
  theme(panel.background = element_rect(fill= "#D1CFD0"), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  scale_y_continuous(breaks = round(c(0,mean(s.a.$pr.po)-0.0065, mean(s.a.$pr.gn)+0.009, 1),2))+
  xlab('iteration')+
  ylab(expression(paste("Mean of  ", Pi("k = K | A"))))

  


plot(s.a.$pr.po, ylim = c(0,1))
abline(h = mean(s.a.$pr.po))
points(s.a.$pr.gn, col = "red", add = T)
abline(h = mean(s.a.$pr.gn), col = "red")
s.a.
