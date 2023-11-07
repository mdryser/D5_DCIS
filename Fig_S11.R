### Create the tree for the tumor

library(knitr)
knitr::opts_chunk$set(echo = FALSE)

library(tidygraph)
library(ggraph)
library(tidyverse)
library(dplyr)
library(igraph)
library(graphlayouts)
library(scatterpie)
library(gplots)
library(Rtsne)
library(scatterplot3d)
library(pracma)
library(diptest)
library(ggExtra)

source('../100_Universal_functions/functs.R')


lambda=1/3 #branching rate; mean branch length (mm) is 1/lambda
G=30 #  max branch level
q=0.5  # Probability of branch termination ("T") vs further branching ("B")

#---------------------------
#---- Example 1 ------------
#---------------------------

set.seed(990)
dum=0 # dummy variable
while(dum<G){
  ductal_tree = ductal.tree.sim(lambda=lambda, G=G, q=q)
  dum=max(vertex_attr(ductal_tree, "level", index = V(ductal_tree)))
}

ductal_tree <-ductal_tree %>%
  activate(edges) %>%
  mutate(colr=ifelse(from<4, 1, 0))

pdf(file = "Tree_example_1.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

ductal_tree %>%
  ggraph(layout = "dendrogram", length=rep(1, gsize(ductal_tree))) +
  geom_node_point() +
  geom_edge_diagonal(aes(colour=colr)) +
  theme(panel.background = element_blank()) 

dev.off()



#---------------------------
#---- Example 2 ------------
#---------------------------

set.seed(9990)
dum=0 # dummy variable
while(dum<G){
  ductal_tree = ductal.tree.sim(lambda=lambda, G=G, q=q)
  dum=max(vertex_attr(ductal_tree, "level", index = V(ductal_tree)))
}

ductal_tree <-ductal_tree %>%
  activate(edges) %>%
  mutate(colr=ifelse(from<4, 1, 0))

pdf(file = "Tree_example_2.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 10) # The height of the plot in inches

ductal_tree %>%
  ggraph(layout = "dendrogram", length=rep(1, gsize(ductal_tree))) +
  geom_node_point() +
  geom_edge_diagonal(aes(colour=colr)) +
  theme(panel.background = element_blank()) 

dev.off()

# ductal_tree %>%
#   ggraph(layout = "dendrogram", length=rep(1, gsize(ductal_tree))) +
#   geom_node_point() +
#   geom_edge_diagonal() +
#   theme(panel.background = element_blank()) 

# ductal_tree %>%
#   ggraph(layout = "dendrogram", length=rep(1, gsize(ductal_tree))) +
#   geom_node_point(size=.8) +
#   geom_edge_link0() + theme_classic()

#------ Visualize the tree as a dendrogram
 ductal_tree %>%
   ggraph(layout = "dendrogram", length=length) +
   geom_node_point() +
   geom_edge_diagonal() +
   geom_node_text(aes(label = number), repel=TRUE)

 ductal_tree %>%
   ggraph(layout = "dendrogram", length=rep(1, gsize(ductal_tree))) +
   geom_node_point() +
   geom_edge_link0() +
   geom_node_text(aes(label = number), repel=TRUE)
 