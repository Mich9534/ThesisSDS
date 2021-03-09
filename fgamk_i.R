# Script for the denominator of the sum V_n(t): (gamma * k)^(n)

fgamk_i <- function(gamma_, i, k){
  
  gk_i <- gamma_ * k 
  
  if(i-1 ==0) {return(gk_i)}
  if(i == 0){ return(1)}
  
  for (j in 1:(i-1)) {
    gk_i <- gk_i * ((gamma_ * k) + j)
    
  }
  return(gk_i)
}


