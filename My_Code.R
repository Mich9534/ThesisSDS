#-----------------------------------------------------------------------------------------------------------------------------------------
#                                                 IMPLEMENTATION GIBBS SAMPLER GENG ET AL.
#-----------------------------------------------------------------------------------------------------------------------------------------
library(doSNOW)
library(ggplot2)
require(parallel)
library(gdata)
library(extraDistr)

#--------------------------------------------------------------------
# RECALL MY FUNCTIONS
#--------------------------------------------------------------------
source('Vn_Miller.R')
source('fv_it.R')
source('fgamk_i.R') # here i = n already implemented into fgamk
source('fk_t.R')  
source('fgenera_data.R')
source('ApproximationVn.R')
source('lVnGn.R')


#-----------------------------------------------------------------------------------------------------------------------------------------
# Let us generate data and see if the function written above works.
#-----------------------------------------------------------------------------------------------------------------------------------------
gamma_ <- 1
n = 100
a = b = 1
data_gen <- fgenera_data(n, 2, 0.3, 0.1, "Balanced")
A <- data_gen$A0
AAA = matrix(0,n,n) ##the upper traiangle for the adjacency matrix
for (i in 1:n){
  for (j in i:n){
    AAA[i,j] = A[i,j]
  }
}
diag(AAA) = 0 ## make it without-selfloop network
diag(A) = 0 ## make it without-selfloop networkk


z0 <- data_gen$z0


start.time <- Sys.time()
logVn.Miller <- log_Vn.M(gamma_, n, n+10)
z1 <- fGibbs_Sampler(100, 9, A, n, a, b, gamma_ = 1, logVn.Miller, gamma.gn = 0.5)  # check if the algorithm works fine
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

start.time <- Sys.time()
z2 <- CDMFM_new_gn(data = A, data1 = AAA, niterations = 100, beta.a = 1, beta.b = 1, GAMMA=1, LAMBDA = 1, initNClusters = 9, gamma.gn = 0.5)
end.time <- Sys.time()
time.taken2 <- end.time - start.time 
time.taken2
# SEMBRA CHE LA MIA IMPLEMENTAZIONE SIA PIÙ RAPIDA 2.4  SEC VS 7.69  SEC con N = 100
# ------------------------------------------------ 2.3  MIN VS 31.51 SEC CON N = 300 -> sembra ci sia parità in questa situazione.. (ho ripetuto più volte l'esperimento)
# ------------------------------------------------ 1.19 MIN VS 2.71  MIN CON N = 500

# VEDIAMO ORA SE LA DIFFERENZA TRA GNEDIN E POISSON PRIOR SI FA SENTIRE
start.time <- Sys.time()
logVn.Miller <- log_Vn.M(gamma_, n, n+10)
z1 <- fGibbs_Sampler(100, 9, A, n, a, b, gamma_ = 1, logVn.Miller, gamma.gn = 0.5)  # check if the algorithm works fine
end.time <- Sys.time()
time.taken <- end.time - start.time

start.time <- Sys.time()
z2 <- fGibbs_Sampler(100, 9, A, n, a, b, gamma_ = 1, logVn.Miller, gamma.gn = 0.5) 
end.time <- Sys.time()
time.taken2 <- end.time - start.time 


table(z)
table(z$Iterates)

#-----------------------------------------------------------------------------------------------------------------------------------------  
# Let's do a parallelization version of the 10 chains run:
#-----------------------------------------------------------------------------------------------------------------------------------------
gamma_ <- 1
n = 100
a = b = 1
p = 0.5
q = 0.1
K = 2
M = 400
K_start = 9

K_tot <- rep(NA, 100)
for(j in 1:100){
  print(paste("dataset N.:",j))
  data_gen <- fgenera_data(n, K, p, q, type_network = "Balanced")    # Let's generate a new dataset with the same specifics
  A <- data_gen$A0
  
  cl <- makeCluster(3)
  registerDoSNOW(cl)
  result <- foreach(1:10) %dopar% fGibbs_Sampler(M, K_start, A, n, a, b, gamma_, logVn.Miller = NA, gamma.gn = 0.5)
  K_chains <- rep(NA, 10)   # Vector where I'll put all my results
  for(i in 1:10){
    nk <- result[[i]][2]
    K_chains[i] <- as.numeric(names(which(table(nk$num_k[-(1:150)]) == max(table(nk$num_k[-(1:150)])))))
  }
  K_chains
  K_tot[j] <- as.numeric(names(which(table(K_chains) == max(table(K_chains)))))
  stopCluster(cl)
}

#--------------------------------------------------------------------------------------------------------------------------------
# Replica of Figure 2 of Geng-Battachya
#--------------------------------------------------------------------------------------------------------------------------------
sample.size <- seq(10, 90, 15)
gamma_ <- 1
a = b = 1
p = 0.5
q = 0.1
K = 2
M = 2000
K_start = 3

prob <- rep(NA, length(sample.size))

## taking the data into the MFM-SBM algorithm
for(s in 1:length(sample.size)){
  n = sample.size[s]
  print(paste("size:", sample.size[s]))
  K_tot <- rep(NA, 100)
  K_chains <- rep(NA, 100)
  for (iter in 1:100) {
    print(paste("dataset N.:",iter))
    data_gen <- fgenera_data(n, K, p, q, type_network = "Balanced")    # Let's generate a new dataset with the same specifics
    A <- data_gen$A0
    
    logVn.Miller <- log_Vn.M(gamma_, n, n+10)   
    fit = fcogibbs_po(M, K_start, A, n, a, b, gamma_, logVn.Miller)
    
    #pr_ch <- as.numeric(which(table(fit$num_k[-c(1:700)]) == max(table(fit$num_k[-c(1:700)]))))
    pr_ch <- mean(fit$num_k[-c(1:100)] == K)
    K_chains[iter] <- pr_ch
  }
  prob[s] <- mean(K_chains)
}
prob
plot(prob~sample.size, type = "l", ylim= c(0,1))
#----------------------------------------------------------------------------------------------------------------------------------------   
# Let's plot now the results:
#-----------------------------------------------------------------------------------------------------------------------------------------
x <- as.numeric(names(table(K_tot)))
counts <- as.vector(table(K_tot))
df <- data.frame(x=x, counts=counts)
    
plt <- ggplot(df) + geom_bar(aes(x=x, y=counts), stat="identity")+scale_x_discrete(name ="N. Groups", 
                                                                                   limits=c("1","2","3","4", "5", "6"))+ylim(0,100)
print(plt)
  
  