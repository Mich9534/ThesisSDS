#------------------------------------------------
# PUT CLUSTER LABELS IN BINARY MATRIX FORM
# needed for point estimation of cluster conf.
#------------------------------------------------

vec2mat <- function(clust_lab){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0,V,H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}
