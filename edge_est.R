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
