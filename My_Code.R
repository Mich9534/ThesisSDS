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

fGibbs_Sampler <- function(M, K_start, A, n, a, b, gamma_, logVn.Miller, gamma.gn){
  #------------------------------------------------------------------------
  # Note: For implementations different from Miller2018 remove logVn.Miller  
  #------------------------------------------------------------------------
  
  library(gdata)
  library(extraDistr)
  num_k <- rep(NA, M)
  z_update <- sample(1:K_start, n, replace = T)
  z_start <- z_update
  Q <- matrix(rep(NA, K_start*K_start), ncol = K_start, nrow = K_start)
  
  for(iter in 1:M){
    print(paste("iteration N.:", iter))
    
    k <- length(unique(z_update))    # number of actual clusters
    num_k[iter] <- as.numeric(k)
    
    for(r in 1:k){
      for(s in r:k){
        a_sum <- sum(A[which(z_update == r), which(z_update == s)])
        a_rs <- ifelse(r == s, 1/2 * a_sum, a_sum)
        n_r <- sum(z_update == r)
        n_s <- sum(z_update == s)
        n_rs <- ifelse(r == s, n_r*(n_r - 1)/2, n_r * n_s)
        
        #print(paste("a_rs + a:", a_rs + a, "  n_rs - a_rs + b:",n_rs - a_rs + b))
        Q[r,s] <- rbeta(1, a_rs + a, n_rs - a_rs + b)
        
      }
    }
   lowerTriangle(Q) <- upperTriangle(Q, byrow = TRUE)  
    
   #print(isSymmetric(Q)) #Checkpoint
   
   
   for(i in 1:n){
     #print(paste("i:",i))  # checkpoint
     pr_tbl <- rep(NA, k)
     for(tbl in 1:k){
       indx <- (1:n)[-i]
       prod <- 1
       for (ind in indx) {
         prod <- prod*(Q[tbl, z_update[ind]])^(A[i, ind]) * (1 - Q[tbl, z_update[ind]])^(1 - A[i, ind])
       }
       #print(paste("prod2:", prod))  # checkpoint
       pr_tbl[tbl] <- (sum(z_update==tbl)+gamma_)* prod
     }
     
     c_i <- length(unique(z_update[-i]))
     
     mAi <- 1
     for(t in 1:c_i){
       j_ct <- which(z_update[-i] == t)
       somma <- sum(A[j_ct, i])
       mAi <- mAi * beta(a, b)^(-1) * beta(somma + a, length(j_ct) - somma + b)
     }
     
     #---------------------------------------------------------------------------------
     # GNEDIN PRIOR
     #---------------------------------------------------------------------------------
     pr_tbl <- c(pr_tbl, exp(l.Vn.Gn(n,c_i+1, gamma.gn)-l.Vn.Gn(n,c_i, gamma.gn))*(gamma_)*mAi)
     
     #---------------------------------------------------------------------------------
     # MILLER IMPLEMENTATION
     #---------------------------------------------------------------------------------
     #pr_tbl <- c(pr_tbl, exp(logVn.Miller[c_i+1]-logVn.Miller[c_i])*(gamma_)*mAi)
     
     #---------------------------------------------------------------------------------
     # OTHER IMPLMENENTATION OF COMPUTATION OF V_n(t). Need to remove Vn.vec 
     #pr_tbl <- c(pr_tbl, exp(flv_n.app(n, c_i+1, gamma_)-flv_n.app(n, c_i, gamma_))*gamma_*mAi) # Using Vn approximated as theorem 5.1    
     #pr_tbl <- c(pr_tbl, (fv_nt(n, c_i+1)/fv_nt(n, c_i)) * gamma_ * mAi) 
     #---------------------------------------------------------------------------------
     
     k.old <- k
     zi.old <- z_update[i]
     z_update[i] <- sample(1:(k+1), 1, replace = FALSE, prob = pr_tbl)
     k <- length(unique(z_update)) 
     
     if(k.old > k){   
       #-----------------------------------------------------
       # devo solamente rimuovere una parte della matrice Q
       #-----------------------------------------------------
       Q <- as.matrix(Q[-zi.old, -zi.old])
       z_update[which(z_update > zi.old)] <- z_update[which(z_update > zi.old)]-1
     }
     
     if(k.old < k){   
       #----------------------------------------------------
       # devo parzialmente aggiornare la matrice Q
       #----------------------------------------------------
       new.prb <- rep(NA, k)
       for (group in 1:(k-1)) {
         a_sum <- sum(A[which(z_update == group), which(z_update == k)])
         n_group <-   sum(z_update == group)
         n_k <-   sum(z_update == k)
         n_gk <-  n_group * n_k
         new.prb[group] <- rbeta(1, a_sum + a, n_gk - a_sum + b)
         
       } 
       
       new.prb[k] <- rbeta(1, a, b)
       Q <- as.matrix(rbind(cbind(Q, new.prb[-k],row.names = NULL), new.prb[],row.names = NULL))
     }
     
     if(max(z_update) > k){ 
       #--------------------------------------------------
       # devo unire le due operazioni precedenti
       #--------------------------------------------------
       new.prb <- rep(NA, max(z_update))
       for (group in (1:max(z_update))[-zi.old]) {
         a_sum <- sum(A[which(z_update == group), which(z_update == max(z_update))])
         n_r <-   sum(z_update == group)
         n_s <-   sum(z_update == max(z_update))
         n_rs <-  n_r * n_s
         new.prb[group] <- rbeta(1, a_sum + a, n_rs - a_sum + b)
       }
       new.prb[max(z_update)] <- rbeta(1, a, b)
       Q.par <- cbind(Q, new.prb[-max(z_update)],row.names = NULL)
       Q <- as.matrix(rbind(Q.par, new.prb[],row.names = NULL))
       Q <- as.matrix(Q[-zi.old, -zi.old])
       z_update[which(z_update > zi.old)] <- z_update[which(z_update > zi.old)]-1
    }
   }
  }
  res <- list(z_update, num_k, z_start)
  names(res) <- c("z_update", "num_k", "z_start")
  return(res)
}

#-----------------------------------------------------------------------------------------------------------------------------------------
# Let us generate data and see if the function written above works.
#-----------------------------------------------------------------------------------------------------------------------------------------
gamma_ <- 1
n = 700
a = b = 1
data_gen <- fgenera_data(n, 2, 0.5, 0.1, "Balanced")
A <- data_gen$A0
z0 <- data_gen$z0
logVn.Miller <- log_Vn.M(gamma_, n, 20)
z <- fGibbs_Sampler(500, 9, A, n, a, b, gamma_ = 1, logVn.Miller, gamma.gn = 0.5)  # check if the algorithm works fine
table(z$num_k)

#-----------------------------------------------------------------------------------------------------------------------------------------  
# Let's do a parallelization version of the 10 chains run:
#-----------------------------------------------------------------------------------------------------------------------------------------
gamma_ <- 1
n = 100
a = b = 1
p = 0.24
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

#-----------------------------------------------------------------------------------------------------------------------------------------   
# Let's plot now the results:
#-----------------------------------------------------------------------------------------------------------------------------------------
x <- as.numeric(names(table(K_tot)))
counts <- as.vector(table(K_tot))
df <- data.frame(x=x, counts=counts)
    
plt <- ggplot(df) + geom_bar(aes(x=x, y=counts), stat="identity")+scale_x_discrete(name ="N. Groups", 
                                                                                   limits=c("1","2","3","4", "5", "6"))+ylim(0,100)
print(plt)
  
