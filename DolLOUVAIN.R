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
diag(A) <- 0
n = dim(A)[1]
net    <- graph.adjacency(A, mode=c("undirected"), weighted=NULL, diag=FALSE)   # transform the adjacency matrix into an igraph object
Louv   <- cluster_louvain(net)$membership                                       # point estimate
K.l    <- length(table(Louv))                                                   # estimated K
VI.l   <- compare(Louv, Z0, method = "vi")
RI.l   <- compare(Louv, Z0, method = "rand")
table(Louv)

# Let us plot the associated block matrix
V = dim(A)[1]
diag(A) <- 0
Y = A
z_0 = Louv
z_1 = z_0[order(z_0)]
Y_1 = Y[order(z_0), order(z_0)]
row_plot_Y_1 <- as.data.frame(as.factor(matrix(z_1,V,1)))
names(row_plot_Y_1) <-"z_1"
rownames(Y_1) <- rownames(row_plot_Y_1)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_Y_1$z_1)
mycolors <- list(z_1 = mycolors)
mycolors[[1]][1] = "#F26D6D"
mycolors[[1]][2] = "#A69B03"
mycolors[[1]][3] = "#04BF7B"
mycolors[[1]][4] = "#00B0F6"
mycolors[[1]][5] = "#E36DF2"
  
Network <- pheatmap(Y_1,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,
                    annotation_row = row_plot_Y_1,annotation_names_row=F,show_rownames=F, show_colnames=F,legend=F,
                    border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)  


# Let us plot the graph
# GGNET
library(ggnetwork)
library(ITNr)
library(intergraph)

net<-ggnetwork(A)
net$vertex.names
g1 = which(Louv  == 1)
g2 = which(Louv  == 2)
g3 = which(Louv  == 3)
g4 = which(Louv  == 4)
g5 = which(Louv  == 5)

col = rep(NA, length(net$vertex.names))
for(i in 1:length(col)){
  if(net$vertex.names[i] %in% g1){
    col[i] = 1
  } else{
    if(net$vertex.names[i] %in% g2){ col[i] =2}else{
      if(net$vertex.names[i] %in% g3){col[i] = 3}else{if(net$vertex.names[i] %in% g4){col[i] = 4} else{col[i] = 5}}
  }
  }
}
net$Group =as.factor(col)
dol.plot = ggplot(net, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(aes(color = Group,size=4)) +
  guides(col = FALSE, size=FALSE)+ 
  theme_blank()



ga = grid.arrange(dol.plot,Network[[4]],nrow=1,widths= c(1.5,1))
g2_PO <- cowplot::ggdraw(ga)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))
