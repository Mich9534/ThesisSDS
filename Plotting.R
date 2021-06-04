# PLOT NETWORK IGRAPH

# Let us plot this with Igraph

library("igraph")


adj <- graph_from_adjacency_matrix(A, mode = "undirected")
node.size= c(10,10,10)

# Node color
n = 120
V(adj)$community <- ifelse(V(adj) <= n/2, 1,2) # for two colors

V(adj)$community <- z0 # for the true confuguration

V(adj)$community <- z$z_update # for the FINAL confuguration

V(adj)$community <- z$z_start # for the starting configuration


rain <- rainbow(14, alpha=0.5)
V(adj)$color <- rain[as.numeric(V(adj)$community)]

E(adj)$color <- apply(as.data.frame(get.edgelist(adj)), 1, 
                      function(x) ifelse(as.numeric(V(adj)$community[x[1]]) == as.numeric(V(adj)$community[x[2]]), 
                                         rain[as.numeric(V(adj)$community[x[1]])], '#E8E8E8'))

E(adj)$width <- apply(as.data.frame(get.edgelist(adj)), 1, 
                      function(x) ifelse(as.numeric(V(adj)$community[x[1]]) == as.numeric(V(adj)$community[x[2]]), 
                                         1, 0.3))

m <- layout_nicely(adj)

plot(adj, vertex.size=6, edge.size = 4, vertex.label=NA, 
     edge.color=E(adj)$color,edge.curved=0.2,
     edge.width=E(adj)$width, layout = m)


# diag(A_0) # no self loops, ok

