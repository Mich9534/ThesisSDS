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