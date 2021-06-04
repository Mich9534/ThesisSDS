rm(list = ls)
#-------------------------------------------------------
# RAND INDEX EXPERIMENT
#-------------------------------------------------------
library(doSNOW)
library(ggplot2)
require(parallel)
library(gdata)
library(extraDistr)
library(fossil)
setwd("/home/michele/Desktop/Small_code_R")
source('Vn_Miller.R')
  source('fv_it.R')
source('fgamk_i.R') # here i = n already implemented into fgamk
source('fk_t.R')  
source('fgenera_data.R')
source('ApproximationVn.R')
source('lVnGn.R')
source('cogibbs_po.R')
source('cogibbs_gn.R')

n = 100
K = 2
p = 0.4
q = 0.1
a = b = gamma_ = 1
gamma.gn = 0.5
par_po = 1
M = 400
N = 100
K_start = 9

set.seed(25)
z_start <- c(sample(1:K_start, size = K_start, replace = FALSE),
             sample(1:K_start, size = n-K_start, replace = TRUE)) # we're sure every cluster is inside furthermore
                                                                  # we want to pass the same starting cluster conf
                                                                  # to all the MCMC to reduce variability of RI

RI.Matrix.gn = matrix(data = rep(NA, M*N), ncol = N)
RI.Matrix.po = matrix(data = rep(NA, M*N), ncol = N)


for(i in 1:N){
  print(paste("iteration number:", i))
  data_gen <- fgenera_data(n, K, p, q, type_network = 'balanced')
  A  = data_gen$A0
  z0 = data_gen$z0
  fit.gn = fcogibbs_gn(M, K_start, A, n, a, b, gamma_, gamma.gn, z0, z_start)
  RI.Matrix.gn[,i] = fit.gn$RI
  
  logVn.Miller = log_Vn.M(gamma_, n, n+10, 1)
  fit.po = fcogibbs_po(M, K_start, A, n, a, b, gamma_, logVn.Miller, z0, z_start)
  RI.Matrix.po[,i] = fit.po$RI
}

meanRI.po <- apply(RI.Matrix.po, 1, mean)
meanRI.gn <- apply(RI.Matrix.gn, 1, mean)

#RI.Matrix
lin.min.po = apply(RI.Matrix.po, 1, min)
lin.max.po = apply(RI.Matrix.po, 1, max)

lin.min.gn = apply(RI.Matrix.gn, 1, min)
lin.max.gn = apply(RI.Matrix.gn, 1, max)

df1 = data.frame(cbind(1:M, lin.min.po, lin.max.po, meanRI.po))
df2 = data.frame(cbind(lin.min.gn, lin.max.gn, meanRI.gn))
df  = cbind(df1, df2)

ggplot(data = df)+
  geom_ribbon(aes(x=V1, ymax=lin.max.po, ymin=lin.min.po), fill="red", alpha=.2)+
  geom_line(aes(x = V1, y = meanRI.po),colour = "red", size = 0.8)+
  geom_ribbon( aes(x=V1, ymax=lin.max.gn, ymin=lin.min.gn), fill="#1FE0CA", alpha=.2)+
  geom_line(aes(x = V1, y = meanRI.gn), colour = "#1FE0CA", size = 0.8)+
  theme(panel.background = element_rect(fill= "white"), panel.border = element_blank(),
                                                                                       panel.grid.major = element_blank(),
                                                                                       panel.grid.minor = element_blank())+
  xlab('Iteration')+
  ylab("Rand Index")+
  geom_segment(aes(x = 300, y = 0.65, xend = 320, yend = 0.65), col = "red")+
  geom_segment(aes(x = 300, y = 0.63, xend = 320, yend = 0.63), col = "#1FE0CA")+
  annotate(geom = "text", x = 340, y = 0.65, label = "Poisson", size = 3)+
  annotate(geom = "text", x = 340, y = 0.63, label = "Gnedin", size = 3)

#-------------------------------------------------------
# Let's try to see if we can do something with dif. par
#-------------------------------------------------------
par_gn <- seq(0.1, 0.9, 0.2)
#par_po <- seq(1,20, 2)
# RI.Matrix.gn1 = matrix(data = rep(NA, M*N), ncol = N)
# RI.Matrix.gn2 = matrix(data = rep(NA, M*N), ncol = N)
# RI.Matrix.gn3 = matrix(data = rep(NA, M*N), ncol = N)
# RI.Matrix.gn4 = matrix(data = rep(NA, M*N), ncol = N)
# RI.Matrix.gn5 = matrix(data = rep(NA, M*N), ncol = N)

RI.Matrix.gn = list(matrix(data = rep(NA, M*N), ncol = N),matrix(data = rep(NA, M*N), ncol = N),matrix(data = rep(NA, M*N), ncol = N),matrix(data = rep(NA, M*N), ncol = N),matrix(data = rep(NA, M*N), ncol = N))

for(i in 1:N){
  print(paste("iteration number:", i))
  data_gen <- fgenera_data(n, K, p, q, type_network = 'balanced')
  A  = data_gen$A0
  z0 = data_gen$z0
  for(ii in 1:length(par_gn)){
    print(paste("par_gn:", par_gn[ii]))
    fit.gn = fcogibbs_gn(M, K_start, A, n, a, b, gamma_, par_gn[ii], z0, z_start)
    RI.Matrix.gn[[ii]][,i] = fit.gn$RI
  }
}

lin.min.gn1 = apply(RI.Matrix.gn[[1]], 1, min)
lin.max.gn1 = apply(RI.Matrix.gn[[1]], 1, max)
lin.min.gn2 = apply(RI.Matrix.gn[[2]], 1, min)
lin.max.gn2 = apply(RI.Matrix.gn[[2]], 1, max)
lin.min.gn3 = apply(RI.Matrix.gn[[3]], 1, min)
lin.max.gn3 = apply(RI.Matrix.gn[[3]], 1, max)
lin.min.gn4 = apply(RI.Matrix.gn[[4]], 1, min)
lin.max.gn4 = apply(RI.Matrix.gn[[4]], 1, max)
lin.min.gn5 = apply(RI.Matrix.gn[[5]], 1, min)
lin.max.gn5 = apply(RI.Matrix.gn[[5]], 1, max)

meanRI.gn1 <- apply(RI.Matrix.gn[[1]], 1, mean)
meanRI.gn2 <- apply(RI.Matrix.gn[[2]], 1, mean)
meanRI.gn3 <- apply(RI.Matrix.gn[[3]], 1, mean)
meanRI.gn4 <- apply(RI.Matrix.gn[[4]], 1, mean)
meanRI.gn5 <- apply(RI.Matrix.gn[[5]], 1, mean)

df.gn = data.frame(cbind(1:M, lin.max.gn1, lin.min.gn1, lin.max.gn2, lin.min.gn2, lin.max.gn3, lin.min.gn3, lin.max.gn4, lin.min.gn4, lin.max.gn5, lin.min.gn5))

ggplot(data = df)+
  geom_ribbon(aes(x=V1, ymax=lin.max.gn1, ymin=lin.min.gn1), fill="#8DEBC3", alpha=.3)+
  geom_ribbon( aes(x=V1, ymax=lin.max.gn2, ymin=lin.min.gn2), fill="#7775EB", alpha=.25)+
  geom_ribbon( aes(x=V1, ymax=lin.max.gn3, ymin=lin.min.gn3), fill="#A7EB6E", alpha=.2)+
  geom_ribbon( aes(x=V1, ymax=lin.max.gn4, ymin=lin.min.gn4), fill="#EB6D57", alpha=.15)+
  geom_ribbon( aes(x=V1, ymax=lin.max.gn5, ymin=lin.min.gn5), fill="#EBD363", alpha=.1)+
  geom_line(aes(x = V1, y = meanRI.gn1),colour = "#8DEBC3", size = 0.8)+
  geom_line(aes(x = V1, y = meanRI.gn2), colour = "#7775EB", size = 0.8)+
  geom_line(aes(x = V1, y = meanRI.gn3), colour = "#A7EB6E", size = 0.8)+
  geom_line(aes(x = V1, y = meanRI.gn4), colour = "#EB6D57", size = 0.8)+
  geom_line(aes(x = V1, y = meanRI.gn5), colour = "#EBD363", size = 0.8)+
  theme(panel.background = element_rect(fill= "white"), panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab('Iteration')+
  ylab("Rand Index")+
  geom_segment(aes(x = 300, y = 0.7,  xend = 320, yend = 0.7 ), col = '#8DEBC3')+
  geom_segment(aes(x = 300, y = 0.68, xend = 320, yend = 0.68), col = '#7775EB')+
  geom_segment(aes(x = 300, y = 0.66, xend = 320, yend = 0.66), col = '#A7EB6E')+
  geom_segment(aes(x = 300, y = 0.64, xend = 320, yend = 0.64), col = '#EB6D57')+
  geom_segment(aes(x = 300, y = 0.62, xend = 320, yend = 0.62), col = '#EBD363')+
  annotate(geom = "text", x = 350, y = 0.7,  label = "gamma = 0.1", size = 3)+
  annotate(geom = "text", x = 350, y = 0.68, label = "gamma = 0.3", size = 3)+
  annotate(geom = "text", x = 350, y = 0.66, label = "gamma = 0.5", size = 3)+
  annotate(geom = "text", x = 350, y = 0.64, label = "gamma = 0.7", size = 3)+
  annotate(geom = "text", x = 350, y = 0.62, label = "gamma = 0.9", size = 3)+
  geom_line(aes(x = V1, y = lin.max.gn1),colour = "#8DEBC3", size = 0.3)+
  geom_line(aes(x = V1, y = lin.min.gn1),colour = "#8DEBC3", size = 0.3)+
  geom_line(aes(x = V1, y = lin.max.gn2),colour = "#7775EB", size = 0.3)+
  geom_line(aes(x = V1, y = lin.min.gn2),colour = "#7775EB", size = 0.3)+
  geom_line(aes(x = V1, y = lin.max.gn3),colour = "#A7EB6E", size = 0.3)+
  geom_line(aes(x = V1, y = lin.min.gn3),colour = "#A7EB6E", size = 0.3)+
  geom_line(aes(x = V1, y = lin.max.gn4),colour = "#EB6D57", size = 0.3)+
  geom_line(aes(x = V1, y = lin.min.gn4),colour = "#EB6D57", size = 0.3)+
  geom_line(aes(x = V1, y = lin.max.gn5),colour = "#EBD363", size = 0.3)+
  geom_line(aes(x = V1, y = lin.min.gn5),colour = "#EBD363", size = 0.3)


