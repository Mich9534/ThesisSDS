# Generate a small dataset with two community, n = 10
'%notin%' <- Negate('%in%')

fgenera_data <- function(n, K, p, q, type_network){
  if(type_network %notin% c('Balanced', 'BALANCED', 'balanced','unbalanced','UNBALANCED', 'Unbalanced')){
    print(paste("Please do insert a valid type_network."))
    print(paste("Available options:"))
    print(paste("Balanced, BALANCED, balanced"))
    print(paste("Unbalanced, UNBALANCED, Unbalanced"))
    return()
  }
  if(type_network %in% c('Balanced', 'BALANCED', 'balanced')){
    'generate the true clustering configuration z_0: K multiple of n'
    
    community_size <- rep(trunc(n/K), K)
    if(trunc(n/K)*K != n){
      community_size[K] <- community_size[K]+(n-sum(community_size))
    }
    community_size
    sum(community_size)
    z0 <- rep(NA, n)
    for(group in 1:K){
      for(i in 1:community_size[group]){
        z0[trunc(n/K)*(group - 1) + i] <- group
      }  
    }
  }
  
  if(type_network %in% c('unbalanced','UNBALANCED', 'Unbalanced')){
    community_size <- rep(trunc(n/K), K)
    
    for(i in 1:(length(community_size)-1)){
      half.i <- floor(community_size[i]/2)
      community_size[i] <- half.i
      community_size[i+1] <- community_size[i+1]+half.i
    }
    community_size[K] <- community_size[K]+(n-sum(community_size))    
    community_size
    z0 <- rep(NA, n)
    for(group in 1:K){
      ifelse(group == 1, start.i <- 1, start.i <- sum(community_size[1:(group-1)])+1)
      for(i in start.i:(community_size[group]+start.i-1)){
        z0[i] <- group
      }  
    }
  }
  
  'construct the matrix Q:'
  Q <- matrix(rep(q, K*K), ncol = K, nrow = K)
  diag(Q) <- rep(p, K)
  
  'construct the adjacency matrix A:'
  A0 <- matrix(rep(0, n*n), n, n)# initialize A_0, the true adjacency matrix
  for(row in 1:n){
    for(col in row:n){
      if(row != col){
        A0[row, col] <- rbinom(1, 1, Q[z0[row], z0[col]])
      }
    }
  }
  lowerTriangle(A0) <- upperTriangle(A0, byrow = TRUE)      
  data_generated <- list(z0, Q, A0)
  names(data_generated) <- c("z0", "Q", "A0")
  return(data_generated)
}


