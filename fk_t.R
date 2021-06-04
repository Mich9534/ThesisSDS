# Script for the numerator of the sum of V_n(t): k_(t)
# (f ahead k_t means we're using a function to compute k_t)

fk_t <- function(t, k){
  
  if(t ==0 || t ==1){
    k_t = max(t*k, 1)
  } else{
    k_t <- k
  
    for(ind in 1:(t-1)){
      k_t <- as.numeric(k_t * (k - ind))
    }
  }
  return(k_t)
}
  
