rm(list = ls())
setwd('/home/michele/Desktop/Small_code_R')
source("cogibbs_gn.R")
source("cogibbs_po.R")
source("fgenera_data.R")
source("fgenera_data.R")
source('Vn_Miller.R')
source('fv_it.R')
source('fgamk_i.R') # here i = n already implemented into fgamk
source('fk_t.R')  
source('fgenera_data.R')
source('ApproximationVn.R')
source('lVnGn.R')

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
#---------------------------------------------------------------------
# !!! CORE PERIPHERY STRUCTURE NETWORK CASE !!!
#---------------------------------------------------------------------

#------------------------------------------------
# PUT CLUSTER LABELS IN BINARY MATRIX FORM  
#------------------------------------------------

vec2mat <- function(clust_lab){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0,V,H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}

#------------------------------------------------
# COMPUTE POSTERIOR CO-CLUSTERING MATRIX  
#------------------------------------------------

pr_cc <- function(z_post){
  # in: posterior sample of assignments (VxN_iter matrix)
  # out: VxV matrix c with elements c[vu]=fraction of iterations in which v and u are in the same cluster
  V <- nrow(z_post)    
  N_iter <- ncol(z_post)
  c <- matrix(0,V,V)
  for (t in 1:N_iter){
    Z <- vec2mat(z_post[,t])
    c <- c + Z%*%t(Z)
  }
  return(c/N_iter)
}

#------------------------------------------------
# COMPUTE MATRIX OF ESTIMATED EDGE PROBABILITIES 
#------------------------------------------------

edge_est <- function(memb,Y,a,b){
  # in: vector of cluster labels (memb), VxV adjancency matrix (Y) and hyperparameters beta priors (a,b)
  # out: matrix of estimated edge probabilities
  z <- dummy(memb)
  H <- ncol(z)
  V <- dim(Y)[1]
  Abs_Freq <- t(z)%*%Y%*%z
  diag(Abs_Freq) <- diag(Abs_Freq)/2
  Tot <- t(z)%*%matrix(1,V,V)%*%z
  diag(Tot) <- (diag(Tot)-table(memb))/2
  Rel_Freq <- (a+Abs_Freq)/(a+b+Tot)
  edge_matr <- z%*%Rel_Freq%*%t(z)
  diag(edge_matr)<-0
  return(edge_matr)
}


#-------------------------------------
# EDGE PROBABILITY MATRIX
#-------------------------------------

#------------------------------------------------
# Generate the block matrix
#------------------------------------------------
V        = 50
K        = 2
p        = 0.5
q        = 0.1
K_start  = 9
a = b    = 1
gamma_   = 1
gamma.gn = 0.5
N_iter   = 15000
burn_in  = 5000

# #k = 5
# block_matrix <- matrix(0.3,5,5)
# block_matrix[1,1] <-0.75
# block_matrix[2,2] <- 0.75
# block_matrix[3,2] <- block_matrix[2,3]<-0.75
# block_matrix[4,4] <- 0.75
# block_matrix[5,4] <- 0.75
# block_matrix[4,5] <- block_matrix[5,4]<-0.75
# block_matrix[2,3] <- block_matrix[3,2]<-0.75
# 
# Q   = block_matrix
# K   = 5
# V   = 200
# nk = V/K
# z_0 = c(rep(1,nk),rep(2,nk),rep(3,nk),rep(4,nk),rep(5,nk))


#K= 2
block_matrix <- matrix(0.3,2,2)
block_matrix[1,1] <-0.4
block_matrix[2,2] <- 0.4
block_matrix[1,2] <- block_matrix[2,1]<-0.1

Q   = block_matrix
K   = 2
V   = 50
nk = V/K
z_0 = c(rep(1,nk),rep(2,nk))

#------------------------------------------------
# Generate the TRUE adjacency matrix Y
#------------------------------------------------

Y <- matrix(0, V, V)# initialize Y, the true adjacency matrix
for(row in 1:V){
  for(col in row:V){
    if(row != col){
      Y[row, col] <- rbinom(1, 1, Q[z_0[row], z_0[col]])
    }
  }
}
lowerTriangle(Y) <- upperTriangle(Y, byrow = TRUE)      

diag(Y) <- 0

#------------------------------------------------
# Plotting the adjacency matrix Y 
#------------------------------------------------

# MATRIX PLOT
row_plot_Y <- as.data.frame(as.factor(matrix(z_0,V,1)))
names(row_plot_Y) <-"z_0"
rownames(Y) <- rownames(row_plot_Y)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_Y$z_0)
mycolors <- list(z_0 = mycolors)
mycolors[[1]][1] = "#F26D6D"
mycolors[[1]][2] = "#00BFC4"
Network1 <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,
                    annotation_row = row_plot_Y,annotation_names_row=F,show_rownames=F, show_colnames=F,legend=F,
                    border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)  


# GGNET

library(ggnetwork)
library(ITNr)
library(intergraph)

#net<-asNetwork(Y) #ggnetwork requires a network object
n<-ggnetwork(Y)
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
  geom_nodes(aes(col = Group, size=4)) +
  guides(col = FALSE, size=FALSE)+ 
  theme_blank()

#K= 2
block_matrix <- matrix(0.3,2,2)
block_matrix[1,1] <-0.4
block_matrix[2,2] <- 0.1
block_matrix[1,2] <- block_matrix[2,1]<-0.4

Q   = block_matrix
K   = 2
V   = 50
nk = V/K
z_0 = c(rep(1,nk),rep(2,nk))

#------------------------------------------------
# Generate the TRUE adjacency matrix Y
#------------------------------------------------

Y <- matrix(0, V, V)# initialize Y, the true adjacency matrix
for(row in 1:V){
  for(col in row:V){
    if(row != col){
      Y[row, col] <- rbinom(1, 1, Q[z_0[row], z_0[col]])
    }
  }
}
lowerTriangle(Y) <- upperTriangle(Y, byrow = TRUE)      

diag(Y) <- 0

#------------------------------------------------
# Plotting the adjacency matrix Y 
#------------------------------------------------

# MATRIX PLOT
row_plot_Y <- as.data.frame(as.factor(matrix(z_0,V,1)))
names(row_plot_Y) <-"z_0"
rownames(Y) <- rownames(row_plot_Y)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[3])
names(mycolors) <- unique(row_plot_Y$z_0)
mycolors <- list(z_0 = mycolors)
mycolors[[1]][1] = "#F26D6D"
mycolors[[1]][2] = "#00BFC4"
Network2 <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,
                     annotation_row = row_plot_Y,annotation_names_row=F,show_rownames=F, show_colnames=F,legend=F,
                     border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)  


# GGNET

library(ggnetwork)
library(ITNr)
library(intergraph)

#net<-asNetwork(Y) #ggnetwork requires a network object
n<-ggnetwork(Y)
n$vertex.names
g1 = which(z_0 == 1)
col = rep(NA, length(n$vertex.names))
for(i in 1:length(col)){
  if(n$vertex.names[i] %in% g1){
    col[i] = 1
  } else{col[i] =2}
}
n$Group =as.factor(col)

dol.plot2 = ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(color = "grey50") +
  geom_nodes(aes(col = Group, size=4)) +
  guides(col = FALSE, size=FALSE)+ 
  theme_blank()



ga2 = grid.arrange(dol.plot,Network1[[4]],dol.plot2,Network2[[4]], nrow=2,widths= c(1.5,1))
g2_PO2 <- cowplot::ggdraw(ga2)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

  
#------------------------------------------------
# Testing the GIBBS sampler with Gnedin prior
#------------------------------------------------

z_start  = c(sample(1:K_start, size = K_start, replace = FALSE),
             sample(1:K_start, size = V-K_start, replace = TRUE))
fit      = fcogibbs_gn(N_iter, K_start, Y, V, a, b, gamma_, gamma.gn, z_0, z_start)
z_post   = fit$z_post


set.seed(1)

# point estimate
c_Z_GN <- pr_cc(z_post[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20) # from package mcclust.ext
memb_Z_GN <- memb_Z_GN_VI$cl


row_plot_GN <- as.data.frame(as.factor(matrix(memb_Z_GN,V,1)))
names(row_plot_GN) <- "memb_Z_GN"
rownames(c_Z_GN) <- rownames(row_plot_GN)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[c(7,3)])
names(mycolors) <- unique(row_plot_GN$memb_Z_GN)
mycolors <- list(memb_Z_GN = mycolors)

Marg <- pheatmap(c_Z_GN,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F,
                 cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,
                 show_colnames=F, legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

g <- grid.arrange(Network[[4]],Marg[[4]],nrow=1,vp=viewport(width=1, height=1))
g2 <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

print(g2)


#------------------------------------------------
# Testing the GIBBS sampler with Poisson prior
#------------------------------------------------

logVn.Miller  = log_Vn.M(gamma_, V, V+10, 1)

fit_po   = fcogibbs_po(N_iter, K_start, Y, V, a, b, gamma_, logVn.Miller, z_0, z_start)

z_post   = fit_po$z_update
  


set.seed(1)

# point estimate
c_Z_GN <- pr_cc(z_post[,(burn_in+1):N_iter])
memb_Z_GN_VI <- minVI(c_Z_GN,method="avg",max.k=20) # from package mcclust.ext
memb_Z_GN <- memb_Z_GN_VI$cl


row_plot_GN <- as.data.frame(as.factor(matrix(memb_Z_GN,V,1)))
names(row_plot_GN) <- "memb_Z_GN"
rownames(c_Z_GN) <- rownames(row_plot_GN)
mycolors <- c(brewer.pal(10,"RdBu")[c(4,7)],brewer.pal(10,"PRGn")[c(7,4)],brewer.pal(9,"YlOrBr")[c(7,3)])
names(mycolors) <- unique(row_plot_GN$memb_Z_GN)
mycolors <- list(memb_Z_GN = mycolors)

Marg <- pheatmap(c_Z_GN,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F,
                 cluster_rows= F,annotation_row = row_plot_GN,annotation_names_row=F,show_rownames=F,
                 show_colnames=F, legend=F,border_color=FALSE,annotation_legend=F,annotation_colors=mycolors)

g <- grid.arrange(Network[[4]],Marg[[4]],nrow=1,vp=viewport(width=1, height=1))
g2 <- cowplot::ggdraw(g)+ theme(plot.background = element_rect(fill=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30)[8]))

print(g2)




