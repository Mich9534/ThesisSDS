rm(list = ls())
setwd("/home/michele/Desktop/Small_code_R")
#-----------------------------------------------------------------
# REAL-WORLD NETWORKS  (RIOLO CANTWELL)
#-----------------------------------------------------------------
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

#-----------------------
# Zachary Karate Club
#-----------------------
data(package = "igraphdata")
data(karate)
Y <- karate$adjacency
z0 = karate$group
z_0 = z0
V   = length(z_0)
n = V
K_start       = 9
a = b         = 1
gamma_        = 1
gamma.gn      = 0.5
lambda        = 2
N_iter        = 2000
burn_in       = 1000
z_start       = c(sample(1:K_start, size = K_start, replace = FALSE),
                  sample(1:K_start, size = n-K_start, replace = TRUE))



diag(Y) <- 0
row_plot_Y <- as.data.frame(as.factor(matrix(z_0,V,1)))
names(row_plot_Y) <-"z_0"
rownames(Y) <- rownames(row_plot_Y)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_Y$z_0)
mycolors <- list(z_0 = mycolors)

Network <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,
                    annotation_row = row_plot_Y,annotation_names_row=F,show_rownames=F, show_colnames=F,legend=F,
                    border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)  

fit.gn <- fcogibbs_gn(N_iter, K_start, Y, V, a, b, gamma_, gamma.gn, z_0, z_start)
hist(fit.gn$num_k[-c(1:burn_in)])

logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
print("entro in fit.po!")
fit.po        = fcogibbs_po(N_iter, K_start, Y, n, a, b, gamma_, logVn.Miller, z_0, z_start)
hist(fit.po$num_k[-c(1:burn_in)], probability = T)
table(fit.po$num_k[-c(1:burn_in)])


'Non ci sono buono risultati in questo caso. Potrebbe essere data dalla scarsa numerosità della rete? Bisognerebbe capire.'

#-----------------------
# American College
#-----------------------

Y <- read.csv(file.choose(), sep=' ', header=F)
Y = Y[,-ncol(Y)]
Network <- pheatmap(Y,cluster_cols = F, cluster_rows= F,
                   annotation_names_row=F,show_rownames=F, show_colnames=F,legend=F,
                    border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)  

fit.gn <- fcogibbs_gn(N_iter, K_start, Y, V, a, b, gamma_, gamma.gn, z_0, z_start)
hist(fit.gn$num_k[-c(1:burn_in)])


logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
print("entro in fit.po!")
fit.po        = fcogibbs_po(N_iter, K_start, Y, n, a, b, gamma_, logVn.Miller, z_0, z_start)

#--------------------
# DOLPHIN DATASET
#--------------------
load("/home/michele/Desktop/Jasa-Article/dolphindata.RData")
Y = A
n = ncol(Y)
K_start = 9
a = b = gamma_ = 1
gamma.gn = 0.5
lambda = 1
z_start       = c(sample(1:K_start, size = K_start, replace = FALSE),
                  sample(1:K_start, size = n-K_start, replace = TRUE))
fit.gn <- fcogibbs_gn(10000, K_start, Y, n, a, b, gamma_, gamma.gn,z0 = NULL,  z_start)
hist(fit.gn$num_k[-c(1:4000)])

logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
fit.po        = fcogibbs_po(10000, K_start, Y, n, a, b, gamma_, logVn.Miller, z0=NULL, z_start)
fit.po$num_k[-c(1:4000)]
hist(fit.po$num_k[-c(1:4000)])
