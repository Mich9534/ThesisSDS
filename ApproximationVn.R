# APPROXIMATION Vn FOR N > 150

flv_n.app <- function(n, t, gamma_){
  gam.t <- gamma_*t
  res <- lfactorial(t)- lfactorial(n) + log(gamma(gam.t)) - log(n^(gam.t-1)) + log(dtpois(t, 1,0))
  
  return(res)
}

#exp(fv_n.app(100000, 2, 1)-fv_n.app(100000, 1, 1))
