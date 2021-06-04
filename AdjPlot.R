rm(list = ls())
setwd('/home/michele/Desktop/Small_code_R')
#-------------------------------------------------------------------------------------
# COMMUNITY ESTIMATION & ADJACENCY MATRIX PLOTTING
# from https://github.com/danieledurante/ESBM/blob/master/Application/application.md
#-------------------------------------------------------------------------------------
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

n             = 200
K             = 5
p             = 0.75
q             = 0.3
DATA          = fgenera_data(n, K, p, q, "unbalanced")
A             = DATA$A0
z0            = DATA$z0
K_start       = 9
a = b         = 1
gamma_        = 1
gamma.gn      = 0.1
lambda        = 7
N_iter        = 2000
burn_in       = 1000
z_start       = c(sample(1:K_start, size = K_start, replace = FALSE),
                    sample(1:K_start, size = n-K_start, replace = TRUE))
#--------------------------
# For the adjacency matrix
#--------------------------
Y   = A
z_0 = z0
V   = n
diag(Y) <- 0
row_plot_Y <- as.data.frame(as.factor(matrix(z_0,V,1)))
names(row_plot_Y) <-"z_0"
rownames(Y) <- rownames(row_plot_Y)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_Y$z_0)
mycolors <- list(z_0 = mycolors)

mycolors[[1]][3] = "#04BFBF"
mycolors[[1]][2] = "#F8766D"
mycolors[[1]][1] = "#EBE59C"

Network <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,
                    annotation_row = row_plot_Y,annotation_names_row=F,show_rownames=F, show_colnames=F,legend=F,
                    border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)  


fit_gn        = fcogibbs_gn(N_iter, K_start, A, n, a, b, gamma_, gamma.gn, z0, z_start)
z_post_gn     = fit_gn$z_post

logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
fit_po        = fcogibbs_po(N_iter, K_start, A, n, a, b, gamma_, logVn.Miller, z0, z_start)
z_post_po     = fit_po$z_post

# ----------------------------------------------------
# EDGE PROBABILITY MATRIX GNEDIN (GN)
# ----------------------------------------------------
#set.seed(1)

# point estimate gn
c_Z_GN <- pr_cc(z_post_gn[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20) # from package mcclust.ext
memb_Z_GN <- memb_Z_GN_VI$cl


row_plot_GN <- as.data.frame(as.factor(matrix(memb_Z_GN,V,1)))
names(row_plot_GN) <- "memb_Z_GN"
rownames(c_Z_GN) <- rownames(row_plot_GN)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[c(3,6,9)])
  names(mycolors) <- unique(row_plot_GN$memb_Z_GN)
mycolors <- list(memb_Z_GN = mycolors)

Marg_GN <- pheatmap(c_Z_GN,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F,
                 cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,
                 show_colnames=F, legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)




# COMBINE THE DIFFERENT FIGURES


g_GN <- grid.arrange(Network[[4]],Marg_GN[[4]],nrow=1,vp=viewport(width=1, height=1))
g2_GN <- cowplot::ggdraw(g_GN)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

print(g2_GN)


# ----------------------------------------------------
# EDGE PROBABILITY MATRIX POISSON (PO)
# ----------------------------------------------------

# point estimate gn
c_Z_PO <- pr_cc(z_post_po[,(burn_in+1):N_iter])
memb_Z_PO_VI <- minVI(c_Z_PO,method="avg",max.k=20) # from package mcclust.ext
memb_Z_PO <- memb_Z_PO_VI$cl


row_plot_PO <- as.data.frame(as.factor(matrix(memb_Z_PO,V,1)))
names(row_plot_PO) <- "memb_Z_PO"
rownames(c_Z_PO) <- rownames(row_plot_PO)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_PO$memb_Z_PO)
mycolors <- list(memb_Z_PO = mycolors)

Marg_PO <- pheatmap(c_Z_PO,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F,
                 cluster_rows= F,annotation_row = row_plot_PO,annotation_names_row=F,show_rownames=F,
                 show_colnames=F, legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)


# COMBINE THE DIFFERENT FIGURES
g_PO <- grid.arrange(Network[[4]],Marg_PO[[4]],nrow=1,vp=viewport(width=1, height=1))
g2_PO <- cowplot::ggdraw(g_PO)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

print(g2_PO)
  
