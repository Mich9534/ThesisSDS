#--------------------------------------------------------------------
# MY IMPLEMENTATION OF COLLAPSED GIBBS SAMPLER USING A GNEDIN PRIOR
#--------------------------------------------------------------------

fcogibbs_gn <- function(M, K_start, A, n, a, b, gamma_, gamma.gn, z_start){
  
  library(gdata)
  library(extraDistr)
  num_k  <- rep(NA, M)
  z_post <- matrix(rep(NA, n*M), ncol = M)
  z_update <- z_start
  Q <- matrix(rep(NA, K_start*K_start), ncol = K_start, nrow = K_start)
  
  for(iter in 1:M){
    #print(paste("iteration N.:", iter))
    
    k <- length(unique(z_update))    # number of actual clusters
    num_k[iter] <- as.numeric(k)
    
    for(r in 1:k){
      for(s in r:k){
        a_sum <- sum(A[which(z_update == r), which(z_update == s)])
        a_rs <- ifelse(r == s, 1/2 * a_sum, a_sum)
        n_r <- sum(z_update == r)
        n_s <- sum(z_update == s)
        n_rs <- ifelse(r == s, n_r*(n_r - 1)/2, n_r * n_s)
        Q[r,s] <- rbeta(1, a_rs + a, n_rs - a_rs + b)
        
      }
    }
    lowerTriangle(Q) <- upperTriangle(Q, byrow = TRUE)  
    
    
    for(i in 1:n){
      
      pr_tbl <- rep(NA, k)
      for(tbl in 1:k){
        
        indx <- (1:n)[-i]
        l.prod <- 0     
        for(ind in c(indx)){
          
          par.sum <- log(Q[tbl, z_update[ind]])*(A[i, ind]) + log(1 - Q[tbl, z_update[ind]])*(1 - A[i, ind])
          
          l.prod <- l.prod+par.sum
          
        }
        pr_tbl[tbl] <- (sum(z_update==tbl)+gamma_) * exp(l.prod)
        #pr_tbl[tbl] <- (sum(z_update==tbl)+gamma_)* exp(loglike(z_update,Q,A,i,n))
        #print(l.pr_tbl[tbl])
      }
      
      
      c_i <- length(unique(z_update[-i]))
      
      mAi <- 1
      for(t in 1:c_i){
        j_ct <- which(z_update[-i] == t)
        somma <- sum(A[j_ct, i])
        mAi <- mAi * beta(a, b)^(-1) * beta(somma + a, length(j_ct) - somma + b)
      }
      
      
      # GNEDIN PRIOR
      
      pr_tbl <- c(pr_tbl, exp(l.Vn.Gn(n,c_i+1, gamma.gn)-l.Vn.Gn(n,c_i, gamma.gn))*(gamma_)*mAi)
      
      
      # MILLER IMPLEMENTATION
      
      #pr_tbl <- c(pr_tbl, exp(logVn.Miller[c_i+1]-logVn.Miller[c_i])*(gamma_)*mAi)
      
      # OTHER IMPLMENENTATION OF COMPUTATION OF V_n(t). Need to remove Vn.vec 
      #pr_tbl <- c(pr_tbl, exp(flv_n.app(n, c_i+1, gamma_)-flv_n.app(n, c_i, gamma_))*gamma_*mAi) # Using Vn approximated as theorem 5.1    
      #pr_tbl <- c(pr_tbl, (fv_nt(n, c_i+1)/fv_nt(n, c_i)) * gamma_ * mAi) 
      
      #print(pr_tbl)
      
      k.old <- k
      zi.old <- z_update[i]
      z_update[i] <- sample(1:(k+1), 1, replace = FALSE, prob = pr_tbl)
      k <- length(unique(z_update)) 
      
      if(k.old > k){   
        Q <- as.matrix(Q[-zi.old, -zi.old])
        z_update[which(z_update > zi.old)] <- z_update[which(z_update > zi.old)]-1
      }
      
      if(k.old < k){  
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
   
    #print(table(z_update))
    z_post[,iter] <- z_update
  }
  
  res <- list(z_post, num_k, z_start)
  names(res) <- c("z_post", "num_k", "z_start")
  return(res)
}
