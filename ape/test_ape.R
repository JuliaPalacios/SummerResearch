
## install CRAN Task View for phylogenetics install.packages('ctv')
## library('ctv') install.views('Phylogenetics')
## update.views('Phylogenetics')

## load ape
library(ape);

## simulate a phylogeny
tree <- rtree(n = 20)
plot(tree, edge.width = 2)

tree
str(tree)

tree <- read.tree(text = "(((A,B),(C,D)),E);")
plot(tree, type = "cladogram", edge.width = 2)

library(apTreeshape)
simulate_kingman(
  
  
)