rm(list = ls())
load("~/Desktop/Jasa-Article/dolphindata.RData")                                           # dolphin data
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


V = dim(A)[1]
diag(A) <- 0
Y = A
z_0 = c(1, 2, 1, 1, 1, 2, 2, 2, 1, 2, 1, 1, 1, 2, 1, 1, 1, 2, 1, 2, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 2, 1, 2, 2, 1, 1, 2, 1)
z_1 = z_0[order(z_0)]
Y_1 = Y[order(z_0), order(z_0)]
row_plot_Y_1 <- as.data.frame(as.factor(matrix(z_1,V,1)))
names(row_plot_Y_1) <-"z_1"
rownames(Y_1) <- rownames(row_plot_Y_1)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_Y_1$z_1)
mycolors <- list(z_1 = mycolors)
mycolors[[1]][1] = "#F8766D"
mycolors[[1]][2] = "#04BFBF"


Network <- pheatmap(Y_1,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,
                    annotation_row = row_plot_Y_1,annotation_names_row=F,show_rownames=F, show_colnames=F,legend=F,
                    border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)  


# GGNET

library(ggnetwork)
library(ITNr)
library(intergraph)

net<-asNetwork(g) #ggnetwork requires a network object
n<-ggnetwork(net)
n$vertex.names
g1 = which(z_0 == 1)
col = rep(NA, length(n$vertex.names))
for(i in 1:length(col)){
  if(n$vertex.names[i] %in% g1){
    col[i] = 1
  } else{col[i] =2}
}
n$Group =as.factor(col)

dol.plot = ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(aes(color = Group,size=4)) +
  guides(col = FALSE, size=FALSE)+ 
  theme_blank()

ga = grid.arrange(dol.plot,Network[[4]],nrow=1,widths= c(1.5,1))
g2_PO <- cowplot::ggdraw(ga)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))



