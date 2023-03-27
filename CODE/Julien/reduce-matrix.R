# Required for graphs
library(igraph)
# Create the matrix
A = matrix(c(0,1,1,
             0,1,0,
             1,1,0),
           nr = 3, byrow = TRUE)
# Create the associated digraph
G = graph_from_adjacency_matrix(A)
# Find the components of G
comp_G = components(G, mode="strong")
# Is the graph strongly connected?
if (comp_G$no==1) {
  writeLines("Matrix is irreducible")
} else{
  # Graph has more than 1 strong component, put the matrix
  # in standard form
  idx = order(comp_G$csize, decreasing = TRUE)
  P = c() # Where we store the rows/columns of P
  for (i in idx) {
    P = c(P, which(comp_G$membership == i))
  }
  # The permutation matrix: obtained by applying P to 
  # columns of the identity matrix
  P = diag(dim(A)[1])[,P]
  # The reduced matrix
  M = t(P)%*%A%*%P
}
