# Function to compute V_n(t). It uses other two functions: fgamk_n & fk_t.

fv_nt <- function(n, t){
  
  if(t != 0){                  # the previous terms of the sum are zero
    k <- t
  } else(k <- 1)
  
  
  
  tol <- 1e-100                # Tollerance we admit -> ask how much we admit
  ans <- 0                     # Answer
  ctrl <- TRUE
  
  while(ctrl == TRUE) {
    num <- fk_t(t, k)
    den <- fgamk_i(gamma_, n, k)
    poi <- dtpois(k, 1, 0)
    
    next_term <- ( num / den ) * poi 
    ans <- ans + next_term
    
    if (abs(next_term)<tol){
      ctrl <- FALSE
    } else{k <- k+1}
  }
  return(ans) 
}




