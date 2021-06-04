rm(list = ls())
load("~/Desktop/Jasa-Article/dolphindata.RData")                                           # dolphin data
setwd('/home/michele/Desktop/Small_code_R')
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


# zero-truncated-poisson-prior with lambda = 1
Z0 = c(1, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 1)
M  = 20000
BI = 10000
a = b = gamma_ = 1
lambda = 1
K_start = 9
n = dim(A)[1]
diag(A) <- 0
logVn.Miller  = log_Vn.M(gamma_, n, n+10, lambda)
z_start <- c(sample(1:K_start, size = K_start, replace = FALSE),
             sample(1:K_start, size = n-K_start, replace = TRUE))
dol.po        = fcogibbs_po(M, K_start, A, n, a, b, gamma_, logVn.Miller, z_start)
c_Z_PO <- pr_cc(dol.po$z_post[,(BI+1):M])

memb_Z_PO_VI <- minVI(c_Z_PO, method="avg",max.k=20) # from package mcclust.ext
memb_Z_PO <- memb_Z_PO_VI$cl
peVI.po <- compare(memb_Z_PO, Z0, method = "vi")
RaIn.po <- compare(memb_Z_PO, Z0, method = "rand")

row_plot_PO <- as.data.frame(as.factor(matrix(memb_Z_PO,n,1)))
names(row_plot_PO) <- "memb_Z_PO"

c_Z_PO = c_Z_PO[order(row_plot_PO), order(row_plot_PO)]
row_plot_PO$memb_Z_PO = row_plot_PO$memb_Z_PO[order(row_plot_PO$memb_Z_PO)]
memb_Z_PO = memb_Z_PO[order(memb_Z_PO)]
rownames(c_Z_PO) <- rownames(row_plot_PO)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[c(3,6,9)])
names(mycolors) <- unique(row_plot_PO$memb_Z_PO)
mycolors <- list(memb_Z_PO = mycolors)
mycolors[[1]][3] = "#5E9FF2"
mycolors[[1]][2] = "#04BF33"
mycolors[[1]][1] = "#F8766D"

Marg_PO <- pheatmap(c_Z_PO,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F,
                    cluster_rows= F,annotation_row = row_plot_PO,annotation_names_row= F,show_rownames=F,
                    show_colnames=F, legend=F, border_color=FALSE, annotation_legend=F, annotation_colors=mycolors)


# GGNET
library(ggnetwork)
library(ITNr)
library(intergraph)

net<-ggnetwork(A)
net$vertex.names
c_Z_PO <- pr_cc(dol.po$z_post[,(BI+1):M])
memb_Z_PO_VI <- minVI(c_Z_PO, method="avg",max.k=20) # from package mcclust.ext
memb_Z_PO <- memb_Z_PO_VI$cl
g1 = which(memb_Z_PO  == 1)
g2 = which(memb_Z_PO  == 2)
col = rep(NA, length(net$vertex.names))
for(i in 1:length(col)){
  if(net$vertex.names[i] %in% g1){
    col[i] = 1
  } else{
    if(net$vertex.names[i] %in% g2){ col[i] =2}else{col[i] = 3}
  }
}
net$Group =as.factor(col)

dol.plot = ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(aes(color = Group,size=4)) +
  guides(col = FALSE, size=FALSE)+ 
  theme_blank()

ga = grid.arrange(dol.plot,Marg_PO[[4]],nrow=1,widths= c(1.5,1))
g2_PO <- cowplot::ggdraw(ga)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

# Gnedin Prior
dol.gn         = fcogibbs_gn(M, K_start, A, n, a, b, gamma_, 0.1,  z_start)

c_Z_GN <- pr_cc(dol.gn$z_post[,(BI+1):M])
memb_Z_GN_VI <- minVI(c_Z_GN, method="avg",max.k=20) # from package mcclust.ext
memb_Z_GN <- memb_Z_GN_VI$cl
peVI.gn <- compare(memb_Z_GN, Z0, method = "vi")
RaIn.gn <- compare(memb_Z_GN, Z0, method = "rand")

row_plot_GN <- as.data.frame(as.factor(matrix(memb_Z_GN,n,1)))

names(row_plot_GN) <- "memb_Z_GN"

c_Z_GN = c_Z_GN[order(row_plot_GN), order(row_plot_GN)]
row_plot_GN$memb_Z_GN = row_plot_GN$memb_Z_GN[order(row_plot_GN$memb_Z_GN)]
memb_Z_GN = memb_Z_GN[order(memb_Z_GN)]
rownames(c_Z_GN) <- rownames(row_plot_GN)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[c(3,6,9)])
names(mycolors) <- unique(row_plot_GN$memb_Z_GN)
mycolors <- list(memb_Z_GN = mycolors)
mycolors[[1]][3] = "#04BFBF"
mycolors[[1]][2] = "#04BF33"
mycolors[[1]][1] = "#F8766D"

Marg_GN <- pheatmap(c_Z_GN,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F,
                    cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row= F,show_rownames=F,
                    show_colnames=F, legend=F, border_color=FALSE, annotation_legend=F, annotation_colors=mycolors)


# GGNET
library(ggnetwork)
library(ITNr)
library(intergraph)

net<-ggnetwork(A)
net$vertex.names
c_Z_GN <- pr_cc(dol.gn$z_post[,(BI+1):M])
memb_Z_GN_VI <- minVI(c_Z_GN, method="avg",max.k=20) # from package mcclust.ext
memb_Z_GN <- memb_Z_GN_VI$cl
g1 = which(memb_Z_GN  == 1)
g2 = which(memb_Z_GN  == 2)
g3 = which(memb_Z_GN  == 3)
col = rep(NA, length(net$vertex.names))
for(i in 1:length(col)){
  if(net$vertex.names[i] %in% g1){
    col[i] = 1
  } else{
    if(net$vertex.names[i] %in% g2){ col[i] =2}else{
      if(net$vertex.names[i] %in% g3){col[i] = 3}else{col[i] = 4}
  }
  }
}
net$Group =as.factor(col)

dol.plot = ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(aes(color = Group,size=4)) +
  guides(col = FALSE, size=FALSE)+ 
  theme_blank()

ga = grid.arrange(dol.plot,Marg_GN[[4]],nrow=1,widths= c(1.5,1))
g2_GN <- cowplot::ggdraw(ga)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

