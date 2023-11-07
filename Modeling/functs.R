# ------------------------------------------------
# ------------------------------------------------
# --- Auxiliary functions
# ------------------------------------------------
# ------------------------------------------------


#---------------------------------------------------
#---------------------------------------------------
#---------------------------------------------------
#---- STOCHASTIC DUCTAL TREE
#---------------------------------------------------
#---------------------------------------------------
#---------------------------------------------------


#---------------------------------------------------
#---- Tree Simulation
#---------------------------------------------------
#-- INPUTS
# lambda = branching rate (mm^--1)
# G = max level depth of the tree
# q = termination probability of branch
#-- OUPUTS
# ductal_tree = the stochastic ductal tree object

ductal.tree.sim <- function(lambda, G, q){
  
  #------  Initialize first three nodes
  node_list <- data.frame(number=c(1,2,3), type=c("B","B","B"), level=c(0,1,1))
  
  #------ Initialize corresponding edges
  edge_list <-data.frame(from=c(1,1), to=c(2,3), length=rexp(2,lambda))
  
  #------ Create the binary tree
  
  for(g in 1:(G-1)){
    
    ind<-which(node_list$type=="B" & node_list$level==g)
    
    for(i in ind){
      
      nn<-dim(node_list)[1]
      
      # Add the two daughter nodes and their T vs B  status
      node_list <- node_list %>%
        add_row(number=nn+1,
                type=ifelse(rbinom(1,1,q)==1,"T", "B"),
                level=g+1) %>%
        add_row(number=nn+2,
                type=ifelse(rbinom(1,1,q)==1,"T", "B"),
                level=g+1)
      
      # Add the two new edges
      edge_list <- edge_list %>%
        add_row(from=node_list$number[i], to=nn+1, length=rexp(1,lambda)) %>%
        add_row(from=node_list$number[i], to=nn+2, length=rexp(1,lambda))
      
    }
    
  }
  
  #------ Transform to a tidygraph object
  ductal_tree <- tbl_graph(nodes = node_list, edges = edge_list, directed = TRUE)
  
  
  #------ Add the depth (in length units) from root
  ductal_tree <- ductal_tree %>% activate(nodes) %>%
    mutate(depth= node_distance_from(
      node_is_root(),
      mode = "out",
      weights = length,algorithm = "automatic"))
  
  
  #----- Find leaves and designate them as terminated
  ductal_tree <- ductal_tree %>% activate(nodes) %>%
    mutate(leaf = node_is_leaf()) %>%
    mutate(type = replace(type, leaf==TRUE, "T"))
  
  return(ductal_tree)
}


# ---------------------------------------------------------
# ---------------------------------------------------------
# SPOT SELECTION ON TREE
# ---------------------------------------------------------
# ---------------------------------------------------------

#--- INPUTS
# ductal_tree: the mammary ductal tree
# spot.sel: whether to do random (R) or slide-based (S) sampling
# spot.sample.n: the number of samples 
# slide.l1, slide.l2, slide.l3: the slide depths for sampling (if spot.sel==S)
# --- OUTPUTS
# - out: vector of the nodes in the ductal tree that will be sampled from during tumor processing

f.spot.selection<-function(ductal_tree, spot.sel, spot.sample.n, slide.l1, slide.l2, slide.l3){
  
  
  if( spot.sample.n> gorder(ductal_tree)) {
    warning('sampling more nodes than available')
    spot.sample.n<- gorder(ductal_tree)
  }
  
  
  if(spot.sel=="R"){
    spot.vec<-sample(gorder(ductal_tree), spot.sample.n)
  } else if (spot.sel=="S"){
    
    node_list_ductal<-ductal_tree %>% 
      activate(nodes) %>%
      as_tibble()
    
    spot.vec.l1<-which(node_list_ductal$level==slide.l1)
    h<-length(spot.vec.l1)
    if (h>round(spot.sample.n/3)){
      spot.vec.l1<-sample(spot.vec.l1, round(spot.sample.n/3))
    }
    
    spot.vec.l2<-which(node_list_ductal$level==slide.l2)
    h<-length(spot.vec.l2)
    if (h>round(spot.sample.n/3)){
      spot.vec.l2<-sample(spot.vec.l2, round(spot.sample.n/3))
    }
    
    spot.vec.l3<-which(node_list_ductal$level==slide.l3)
    h<-length(spot.vec.l3)
    if (h>round(spot.sample.n/3)){
      spot.vec.l3<-sample(spot.vec.l3, round(spot.sample.n/3))
    }
    
    spot.vec<-c(spot.vec.l1, spot.vec.l2, spot.vec.l3)
    
  } else {
    warning('WARNING: invalid spot selection choice')
  }
  
  
  return(spot.vec)
  
}





#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------
#---- TUMOR GROWTH
#-----------------------------------------------------------
#-----------------------------------------------------------
#-----------------------------------------------------------

# ---------------------------------------------------------
# ---------------------------------------------------------
# INITIAL TUMOR EXPANSION
# ---------------------------------------------------------
# ---------------------------------------------------------

#--- INPUTS
# N.SC: Number of stem cells to create
# N.mut: Expected number of mutations
# N.gen: Number of cell divisions for the mutation burst
#--- OUTPUTS
# mut_track: array of dimension N.SC x m, where m are the number of mutations



BigBang.expansion <- function(N.SC, N.mut, N.gen){
  
  if(2^N.gen>N.SC) {warning('creating more stem cells than needed')}
  
  # The Poisson rate rate needed to get the expected number of mutation
  lambda=N.mut/(2^(N.gen+1)-1)
  
  mut_track<-matrix(0, 1,N.mut*3) # array to track the mutations in SCs
  
  # WATCH OUT: we make it a minimum of one mutation in the founder cell
  mut<-1+rpois(1, lambda) # add  mutations to founder cell
  if(mut>0){
    mut_track[1,1:mut]<-1
  }
  
  pos=mut # the position of the most recently added mutation
  
  # Loop over all generations
  for(k in 1:N.gen){
    
    # Duplicate each cell
    mut_track<- rbind(mut_track, mut_track) 
    
    # Number of mutations for each new daughter cell
    mut<-rpois(2^k, lambda)
    
    # Add the new mutations after 'pos'
    for(j in 1:2^k){
      
      if(mut[j]>0){
        
        mut_track[j, (pos+1):(pos+mut[j]) ]<-1
        
        pos=pos+mut[j]
        
      }
      
    }
    
  }
  
  #-- Remove mutation slots that remained empty
  ind<-which(colSums(mut_track)>0)
  mut_track <- mut_track[, ind]
  
  #-- At the end of the N.gen generations of exponential growth, expand
  #-- stochastically to the size of the TEB stem cell pool
  if((N.SC-2^N.gen)>0){
    for(j in 1: (N.SC-2^N.gen)){
      mut_track<- rbind(mut_track,
                        mut_track[sample(dim(mut_track)[1],1),])
    }
  }
  
  
  mut_track<-mut_track[,order(colSums(mut_track),decreasing=TRUE)]
  
  return(mut_track)
}


# ---------------------------------------------------------
# ---------------------------------------------------------
# DCIS GROWTH ON TREE
# ---------------------------------------------------------
# ---------------------------------------------------------

#--- INPUTS
# ductal_tree: the mammary ductal tree
# N: number of caner stem cells in each branch
# start_level: depth (in generations) of the founding node
# type: which DCIS growth dynamics to be used, see below for details
# --- OUTPUTS
# cell_track: array 'no nodes' x 'no stem cells / TEB'
#             columns 1:N/2 are border cells
#             columns N/2+1: N are tip cell
# dcis_tree: augmented ductal tree with tumor annotations
# --- DETAILS for the type:
# --- Cancer SC in TEB only:
# TEB.BT.mix: border-tip separation, high diffusion
# TEB.BT.nomix: border-tip separation, low diffusion 
# TEB.noBT.mix: no border-tip separation, high diffusion
# --- All cells are stem cells
# Branch.mix: exponential growth, high diffusion in branch
# Branch.surface: exponential growth at front only
# Branch.volume: exponential growth in entire volume



dcis.growth <- function(ductal_tree, N, start_level, type, cell.diam, n.duct){
  
  
  #------- Initialize cell tracking matrix
  M=gorder(ductal_tree) # number of nodes in the tree
  cell_track<-matrix(0, M,N) # array to track the cells
  
  #------- Extract the node list
  # Add STATUS variable for tracking of tumor growth
  # 0: tumor has not passed yet
  # 1: tumor has arrived but not yet moved beyond
  # 2: tumor has passed through
  
  node_list<-ductal_tree %>% 
    activate(nodes) %>%
    mutate(status=0) %>%
    as_tibble()
  
  #------- Extract edge list
  # Add a DIRECTION variable for direction of tumor growth 
  # 1: moves from root to leaves
  # -1: moves from leaves to root
  edge_list <- ductal_tree %>%
    activate(edges) %>%
    mutate(direction=1) %>%
    as_tibble()
  
  # ---- Sample the starting node
  if(start_level=="R"){
    strt<-sample(node_list$number,1)  # random starting location
  } else {
    strt<-sample(which(node_list$level==start_level),1) # at level=start_level
  }
  # ---- If it is a B node: downward growth is needed; if D node, no further action
  node_list$status[strt]<-ifelse(node_list$type[strt]=="B", 1, 2)
  
  # ----- Create the celllular composition of the starting node
  # ----- Give each cell a unique ID for further tracking
  # ----- 1:N/2 are border cells
  # ----- N/2+1:N are tip cells
  cell_track[strt, ]<-c(1:N)
  

  
  #---------------------------------------------------
  # ----- Move up to  root (direction of growth: -1)
  #---------------------------------------------------
  
  #----- Find shortest path from start to root (node 1)
  path_p<-shortest_paths(graph=ductal_tree,
                         from=1,
                         to=strt)$vpath[[1]]
  path_p<- rev(path_p)
  
  #----- All edges to the root are traveled in direction -1
  #----- Not needed if we start at the root
  if(length(path_p)>1){ 
    
    for(k in 1:(length(path_p)-1)){
      ind<-which(edge_list$from==path_p[k+1] & edge_list$to==path_p[k])
      edge_list$direction[ind]<--1
    }
    #----- First node up inherits the same stem cell pool
    cell_track[path_p[2],]<-cell_track[path_p[1],]
  }
  
  #----- Travel along rest of the path, until root
  #----- Only needed if ≥1 node between start and root
  if( length(path_p)>2){ # This part only if at least one node between start and root
    
    for(k in 2: (length(path_p)-1)){
      
      P<-path_p[k+1] #  parent node
      D<-neighbors(ductal_tree, path_p[k], mode="out") # both children
      D<-D[which(!(D %in% path_p))] # the so far unvisited child
      

      
      
      # Branching dynamics
      if(type=="TEB.BT.mix"){
        cell_track<- branch.TEB.mix(cell_track, path_p[k], P, D)
      } else if (type=="TEB.BT.nomix"){
        cell_track<- branch.TEB.nomix(cell_track, path_p[k], P, D)
      } else if (type=="TEB.noBT.mix"){
        cell_track<- branch.TEB.mix(cell_track, path_p[k], P, D)
      } else if (type=="Branch.mix"){
        cell_track<- branch.duct.mix(cell_track, path_p[k], P, D)
      } else if (type=="Branch.surface"){
        # Figure out the length of the two ducts
        ind<-which(edge_list$from==P & edge_list$to==path_p[k])
        le1<-round(edge_list$length[ind]/cell.diam)
        ind<-which(edge_list$from==path_p[k] & edge_list$to==D)
        le2<-round(edge_list$length[ind]/cell.diam)
        rm(ind)
        cell_track<- branch.duct.surface(cell_track, path_p[k], P, D,le1,le2, n.duct)
      } else if (type=="Branch.volume"){
        # Figure out the length of the two ducts
        ind<-which(edge_list$from==P & edge_list$to==path_p[k])
        le1<-round(edge_list$length[ind]/cell.diam)
        ind<-which(edge_list$from==path_p[k] & edge_list$to==D)
        le2<-round(edge_list$length[ind]/cell.diam)
        rm(ind)
        cell_track<- branch.duct.volume(cell_track, path_p[k], P, D,le1,le2, n.duct)
      } else {
        warning('Invalid branching dynamics')
      }
      
      #----- This node will not need to be revisited
      node_list$status[path_p[k]]<-2
      #----  The daughter is either a status 2 if "T", or a 1 if "B"
      node_list$status[D]<-ifelse(node_list$type[D]=="T",2,1)
      
    }
  }
  
  
  #---------------------------------------------------------
  # ----- Navigate the root
  #---------------------------------------------------------
  #-----  Only if we do not start at the root
  if(length(path_p)>1){ 
    
    D<-neighbors(ductal_tree, 1, mode="out")
    D<-D[which(!(D %in% path_p))]
    cell_track[D,]<-cell_track[1,]
    node_list$status[D]<-ifelse(node_list$type[D]=="T",2,1)
    
    #------ Root is now finished
    node_list$status[1]<-2 
  }
  
  #---------------------------------------------------------
  # ----- Move down to the leaves (direction of growth: +1)
  #---------------------------------------------------------
  
  #----- Grow tumor downward now (Loop over status=1 nodes)
  #----- At this point all nodes are "equal", no special cases left
  while(any(node_list$status==1)){
    
    #-- pick one that is unfinished
    pos<-which(node_list$status==1)[1]
    
    #-- Its two children
    D1<-neighbors(ductal_tree, pos, mode="out")[1]
    D2<-neighbors(ductal_tree, pos, mode="out")[2]
    

    
    
    # Branching dynamics
    if(type=="TEB.BT.mix"){
      cell_track<- branch.TEB.mix(cell_track, pos, D1, D2)
    } else if (type=="TEB.BT.nomix"){
      cell_track<- branch.TEB.nomix(cell_track, pos, D1, D2)
    } else if (type=="TEB.noBT.mix"){
      cell_track<- branch.TEB.mix(cell_track, pos, D1, D2)
    } else if (type=="Branch.mix"){
      cell_track<- branch.duct.mix(cell_track, pos, D1, D2)
    } else if (type=="Branch.surface"){
      # Figure out the length of the two ducts
      ind<-which(edge_list$from==pos & edge_list$to==D1)
      le1<-round(edge_list$length[ind]/cell.diam)
      ind<-which(edge_list$from==pos & edge_list$to==D2)
      le2<-round(edge_list$length[ind]/cell.diam)
      rm(ind)
      cell_track<- branch.duct.surface(cell_track, pos, D1, D2, le1,le2, n.duct)
    } else if (type=="Branch.volume"){
      # Figure out the length of the two ducts
      ind<-which(edge_list$from==pos & edge_list$to==D1)
      le1<-round(edge_list$length[ind]/cell.diam)
      ind<-which(edge_list$from==pos & edge_list$to==D2)
      le2<-round(edge_list$length[ind]/cell.diam)
      rm(ind)
      cell_track<- branch.duct.volume(cell_track, pos, D1, D2, le1, le2,n.duct)
    } else {
      warning('Invalid branching dynamics')
    }
    
    
    #--- Update the status for current node and its children
    node_list$status[pos]<-2
    node_list$status[D1]<-ifelse(node_list$type[D1]=="B", 1, 2)
    node_list$status[D2]<-ifelse(node_list$type[D2]=="B", 1, 2)
    
  }
  
  # Assemble the tree
  dcis_tree <- tbl_graph(nodes = node_list, edges = edge_list, directed = TRUE)
  
  return(list(dcis_tree, cell_track, strt))
  
}


# ---------------------------------------------------------
# ---------------------------------------------------------
# BRANCHING DYNAMICS
# ---------------------------------------------------------

# ---------------------------------------------------------
#---- INPUTS
#   cell_track: matrix of M TEBs (nodes) by N stem cells
#   pos: current node
#   D1: daughter TEB 1
#   D2: daughter TEB 2
#---- OUTPUTS:
#   cell_track: updated matrix

# ---------------------------------------------------------
# -- branch.TEB.mix
# ---------------------------------------------------------
#-  WHAT: Border-tip dynamics as proposed by Scheele et al.
#-  MECHANIM: Cancer stem cells in TEB; Border-tip separation;
#-  High diffusion leads to random distribution during splitting

branch.TEB.mix <- function(cell_track, pos, D1, D2){
  
  # number of cells in TEB
  N <-dim(cell_track)[2]
  # shuffle the current TEB
  cell.now<-sample(cell_track[pos,])
  # allocate  half to each of the two daughters
  cell.D1<-cell.now[1:(N/2)] 
  cell.D2<-cell.now[(N/2+1):N]
  
  # Duplicate the N/2 cells transferred to the daughters
  cell_track[D1,]<-sample(c(cell.D1, cell.D1))
  cell_track[D2,]<-sample(c(cell.D2, cell.D2))
  
  return(cell_track)
  
}



# ---------------------------------------------------------
# -- branch.TEB.nomix
# ---------------------------------------------------------
#-  WHAT: The border tip dynamics with low diffusion in the TEB;
#-  the non-mixing case in Scheele et al. (see their appendix)
#-  MECHANIM: Cancer stem cells in TEB; Border-tip separation;
#-  Low diffusion leads to spatial segregation of cells in TEB during branching

branch.TEB.nomix <- function(cell_track, pos, D1, D2){
  
  # number of cells in TEB
  N<-dim(cell_track)[2]
  
  # determine the cutting point of the torus
  cell.now<-circshift(cell_track[pos,], sample(0:N-1,1))
  
  # two half-torus' with empty spots to be filled
  cell.D1<-c(cell.now[1:(N/2)], rep(0, N/2))
  cell.D2<-c(cell.now[(N/2+1):N], rep(0,N/2))
  
  # create the k-vector from 0 to N/2-1
  vec<-c(0:(N/2-1))
  
  # use n_{2k+1}=n_{k+1}, for k=0 to N/2-1 
  cell_track[D1,2*vec+1]<-cell.D1[vec+1]
  cell_track[D1,2*vec+2]<-cell.D1[vec+1]
  
  cell_track[D2,2*vec+1]<-cell.D2[vec+1]
  cell_track[D2,2*vec+2]<-cell.D2[vec+1]  
  
  return(cell_track)
  
}



# ---------------------------------------------------------
# -- branch.duct.mix
# ---------------------------------------------------------
#-  WHAT: No cancer stem cells (all cells stem cells)
#-  MECHANIM: Growth is duct by duct, with seeding to the daughters;
#-  high diffusion (well-mixed) within each duct


branch.duct.mix <- function(cell_track, pos, D1, D2){
  
  N <- dim(cell_track)[2]
  # Sample with replacement before distributing half-half
  cell.now<-sample(cell_track[pos,], replace = TRUE)
  
  cell_track[D1,]<-c(cell.now[1:(N/2)], cell.now[1:(N/2)])
  cell_track[D2,]<-c(cell.now[(N/2+1):N], cell.now[(N/2+1):N])
  
  return(cell_track)
  
}


# ---------------------------------------------------------
# -- branch.duct.surface
# ---------------------------------------------------------
#-  WHAT: No cancer stem cells (all cells stem cells)
#-  MECHANIM: Growth is duct by duct, with seeding to the daughters;
#-  low-diffusion within each duct, growth at  expanding front/surface only


branch.duct.surface<- function(cell_track, pos, D1, D2, le1, le2, n.duct){
  
  N <- dim(cell_track)[2]
  
  # Split half-half to the beginning of the duct
  cell.now<-sample(cell_track[pos,])
  cell.D1<-c(cell.now[1:(N/2)], cell.now[1:(N/2)])
  cell.D2<-c(cell.now[(N/2+1):N], cell.now[(N/2+1):N])
  
  # sample with replacement to propagate the surface
  t1<-max(1,round(le1*n.duct/N))
  for(k in 1:t1){
    cell.D1<-sample(cell.D1, replace=TRUE)
  }
  
  t2<-max(1,round(le2*n.duct/N))
  for(k in 1:t2){
    cell.D2<-sample(cell.D2, replace=TRUE)
  }
  
  
  cell_track[D1,]<-cell.D1
  cell_track[D2,]<-cell.D2
  
  return(cell_track)
  
  
}


# ---------------------------------------------------------
# -- branch.duct.volume
# ---------------------------------------------------------
#-  WHAT: No cancer stem cells (all cells stem cells)
#-  MECHANIM: Growth is duct by duct, with seeding to the daughters;
#-  low-diffusion within each duct, growth in the entire duct volume


branch.duct.volume <- function(cell_track, pos, D1, D2, le1, le2, n.duct){
  
  N<-dim(cell_track)[2]
  
  # Split half-half to the beginning of the duct
  cell.now<-sample(cell_track[pos,])
  cell.D1<-c(cell.now[1:(N/2)], cell.now[1:(N/2)])
  cell.D2<-c(cell.now[(N/2+1):N], cell.now[(N/2+1):N])
  
  # Update le1 and le2 to account for the volume growth behind the surface
  t1<-max(1,round(log2(le1*n.duct/N)))
  t2<-max(1, round(log2(le2*n.duct/N)))
  
  # sample with replacement to propagate the surface
  for(k in 1:t1){
    cell.D1<-sample(cell.D1, replace=TRUE)
  }
  
  for(k in 1:t2){
    cell.D2<-sample(cell.D2, replace=TRUE)
  }
  
  
  cell_track[D1,]<-cell.D1
  cell_track[D2,]<-cell.D2
  
  return(cell_track)
  
}




# ---------------------------------------------------------
# ---------------------------------------------------------
# CREATE MAF MATRIX
# ---------------------------------------------------------
# ---------------------------------------------------------
#---- INPUTS
#- cell_track: matrix (nodes by stem cells) of which cells in which node
#- mut_track: matrix (stem cells by mutations) of which mutations in which cell
#- type: type  of branching used
#- patch_size: 2^nTA, where nTA=number of transit-amplifications
#- n.duct: number of cells in the duct cross-section
#---- OUTPUT
#- matrix (no mutations x no nodes) with the MAFs

MAF.conv <- function(cell_track, mut_track, type, patch_size, n.duct){
  
  # effective number of cells that create the duct
  N.eff = ifelse(type %in% c("TEB.BT.mix","TEB.BT.nomix"), dim(cell_track)[2]/2, dim(cell_track)[2])
  
  # number of nodes
  n.nodes<-dim(cell_track)[1]
  
  # number of mutations
  n.mut<-dim(mut_track)[2]
  
  # output matrix
  out<-array(0, c( n.mut, n.nodes))
  
  # # The OLD loop
  # for(k in 1:n.mut){ # loop over the mutations
  #   # find cells that contain the mutation
  #   cells.w.mut<-which(mut_track[,k]==1)
  #   # find which nodes contain the mutation
  #   for(m in 1: n.nodes){ # loop over the nodes
  #     # new version with patch_size
  #     dum<-sample(N.eff, max(1, ceiling(N.eff/patch_size)), replace = FALSE)
  #     dum<-rep(dum, patch_size)
  #     out[k,m]<-sum(cell_track[m,dum]%in%cells.w.mut)/N.eff
  #   }
  # }
  
  # number of cells that need to be sampled to make up the duct
  sadu<-ceiling(n.duct/patch_size)
  # The NEW LOOP
  for(m in 1: n.nodes){ # loop over the nodes
    # sample sadu cells from the N.eff cells with replacement
    dum<-sample(N.eff, sadu, replace = TRUE)
    # scale back up to n.duct cells (rounding errors...)
    dum<-rep(dum, patch_size)
    # pick only n.duct of those
    dum<-sample(dum, n.duct, replace=FALSE)
    
    for(k in 1:n.mut){ # loop over the mutations
      # find cells that contain the mutation
      cells.w.mut<-which(mut_track[,k]==1)
      # calculate their frequency in the sample of n.duct cells
      out[k,m]<-sum(cell_track[m,dum]%in%cells.w.mut)/length(dum)
    }
  }
  
  return(out)
  
}



# ---------------------------------------------------------
# ---------------------------------------------------------
# MUTATION DIVERGENCE 1
# ---------------------------------------------------------
# ---------------------------------------------------------
#---- INPUTS
#- mat: matrix where columns are the vectors 
#- ind: among which columns to calculate the divergence
#---- OUTPUTS
#- div: average of the pairwise L1-distance (normalized by length of vector)

f.div.single <- function(mat, ind){
  
  n<-length(ind)
  div<-0
  
  for(q1 in 1:(n-1)){
    for(q2 in (q1+1):n){
      div<-div+sum(abs(mat[,ind[q1]]-mat[,ind[q2]]))/dim(mat)[1]
    }
  }
  
  div<- div/(n*(n-1)/2)
  return(div)
  
}

# ---------------------------------------------------------
# ---------------------------------------------------------
# MUTATION DIVERGENCE 2
# ---------------------------------------------------------
# ---------------------------------------------------------
#---- INPUTS
#- mat: matrix where columns are the vectors 
#- ind1, ind2: among which columns to calculate the divergence
#---- OUTPUTS
#- div: average of the pairwise L1-distance (normalized by length of vector)
f.div.pair <- function(mat, ind1, ind2){
  
  n1<-length(ind1)
  n2<-length(ind2)
  
  div<-0
  
  for(q1 in 1:n1){
    for(q2 in 1:n2){
      div<-div+sum(abs(mat[,ind1[q1]]-mat[,ind2[q2]]))/dim(mat)[1]
    }
  }
  
  div<- div/(n1*n2)
  return(div)
  
}




# ---------------------------------------------------------
# ---------------------------------------------------------
# Genetic-spatial correlation
# ---------------------------------------------------------
# ---------------------------------------------------------
#---- INPUTS
#- mat: matrix where columns are the vectors containing genetic info for each spot
#- coord: vector containing the depth of the nodes
#---- OUTPUTS
#- rho: pearson's r for the genetic vs spatial distance of the vectors

f.spat.gen.corr <- function(mat, coord){
  
  # number of nodes
  n<-length(coord)
  # spatial distance matrix 
  dist.spat<-  as.matrix(dist(coord))
  # prepare genetic distance matrix
  dist.gen<-  NA*dist.spat
  
  # calculate pairwise L1-distance between the spots
  for(q1 in 1:(n-1)){
    for(q2 in (q1+1):n){
      dist.gen[q1,q2]<-sum(abs(mat[,q1]-mat[,q2]))
    }
  }
  
  # calculate Pearson's R
  ind<-upper.tri(dist.gen, diag = FALSE)
  x<-as.vector(dist.gen[ind])
  y<-as.vector(dist.spat[ind])
  rho<-cor(x,y, method="pearson")
  
  rho<-ifelse(is.na(rho), 1, rho)
  
  return(rho)
  
  
}



# ---------------------------------------------------------
# ---------------------------------------------------------
# Expansion index function
# ---------------------------------------------------------
# ---------------------------------------------------------
#---- INPUTS
#- mat: dataframe with two columns for mutation diameter and fraction of spots
#- x
#---- OUTPUTS
#- out: the function value at x

f.EI <- function(mat,x){
  
  ind<-which(mat$spots<=x)
  
  if(length(ind)){
    out<-max(mat$diam[ind])
  } else {
    out<-0
  }
  
  return(out)
  
}

f.EI<-Vectorize(f.EI, vectorize.args = "x")


# ---------------------------------------------------------
# ---------------------------------------------------------
# Expansion index calculation
# ---------------------------------------------------------
# ---------------------------------------------------------
#---- INPUTS
#- MUT.matrix: the mutation matrix (mutations X nodes) of selected nodes
#- node_list_sel: the nodel_list entries from the selected nodes
#---- OUTPUTS
#- mut.stat: dataframe with two columns: normalized diameter and fraction spots covered 

f.expansion.index <- function(MUT.matrix, node_list_sel){
  
  # To record the diameter and spot fraction statistics
  mut.stat<-data.frame(diam=rep(NA, dim(MUT.matrix)[1]),
                       spots=rep(NA, dim(MUT.matrix)[1]))
  
  # Loop over all mutations
  for(j in 1: dim(MUT.matrix)[1]){
    # find spots that contain the mutation
    ind<-which(MUT.matrix[j,]==1)
    # if there are at least two spots...
    if(length(ind)>1){
      # record their respective depth
      depth<-node_list_sel$depth[ind]
      # diameter these spots span divided by the total sampled tumor diameter
      mut.stat$diam[j] <- (max(depth)-min(depth))/max(node_list_sel$depth)
      # fraction of spots that contain the mutaiton among all sampled spots
      mut.stat$spots[j]<-length(ind)/dim(MUT.matrix)[2]
    }
    
  }
  
  # Remove mutations that are in less than two spots
  mut.stat<-mut.stat[which(!is.na(mut.stat$diam)),]
  
  # The actual expansion index
  EI<-integrate(f.EI, lower=0, upper=1, mat=mut.stat)$value
  out<-list(mut.stat,EI)

  return(out)
  
}



# ---------------------------------------------------------
# ---------------------------------------------------------
# Ising energy calculation
# ---------------------------------------------------------
# ---------------------------------------------------------
#---- INPUTS
#- The mutation matrix (already clustered): mutations X nodes
#---- OUTPUTS
#- data frame with one column (energy) that contains the ising energy of each present mutation

f.ising <- function(MUT.matrix){
  
  # Only pick mutations that are present in the sampled tumor
  ind<-which(rowSums(ENERGY.matrix)>0)
  energy.vec <- rep(0,length(ind))
  
  # Calculate the Ising energy for the mutation
  for(k in 1:length(ind)){
    energy.vec[k]<-sum(abs(diff(ENERGY.matrix[ind[k],])))
  }
  
  # Normalize it by the maximum possible number of flips
  energy.vec<-data.frame(energy=energy.vec/(dim(ENERGY.matrix)[2]-1)) 
  
}


#### MULTIPLOT

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}





# ---------------------------------------------------------
# ---------------------------------------------------------
# Prepare the experimental data for further analysis
# ---------------------------------------------------------
# ---------------------------------------------------------
#- INPUT:
#-- tumor.exclude: vector of integers with the labels of tumors to exclude
#- OUTPUT:
#-- Curated data sheets for M1-M6


f.tumor.prep<-function(tumor.exclude){
  
  path_dat_sum<-"/Users/mdr30/Box/0_D5_Paper/6_Modeling/2_Data_summaries"
  
  #-- M1: total number of mutations
  #-- M2: fraction of public mutations
  data_M1_M2<-read.csv(file=paste(path_dat_sum, "/data_M1_M2.csv", sep=""))
  
  data_M1_M2 <- data_M1_M2 %>%
    mutate(n.mutations=mut.total,
           n.public=mut.public,
           n.private=mut.private,
           type="data") %>%
    filter(! (tumor %in% tumor.exclude)) %>%
    select(n.mutations, n.public, n.private, type)
  
  
  #-- M3: Pearson's R for genetic vs spatial spot distance
  data_M3<-read.csv(file=paste(path_dat_sum,"/data_M3.csv", sep=""))
  data_M3 <- data_M3 %>%
    mutate(spatial.corr=Mean,
           type="data") %>%
    filter(! (Tumor %in% tumor.exclude) )%>%
    select(spatial.corr, type)
  
  #-- M4: Expansion index
  data_M4<-read.csv(file=paste(path_dat_sum, "/data_M4.csv", sep=""))
  data_M4 <- data_M4 %>%
    mutate(EI=Mean,
           type="data")%>%
    filter(! (Tumor %in% tumor.exclude) )%>%
    select(EI, type)
  
  #-- M5: Ising Energy median and IQR
  data_M5<-read.csv(file=paste(path_dat_sum, "/data_M5.csv", sep=""), header = TRUE)
  data_M5 <- data_M5 %>%
    select(-X) %>%
    mutate(energy.median= Median) %>%
    mutate(energy.IQR=IQR) %>%
    mutate(type= "data") %>%
    filter(! (Tumor %in% tumor.exclude)) %>%
    select(energy.median, energy.IQR, type)
  
  #-- M6 MAF
  dat.temp<-read.csv(file=paste(path_dat_sum, "/data_M6.csv", sep=""))
  tum<-unique(dat.temp$tumor)
  tum<-tum[!tum %in% tumor.exclude]
  data_M6<-data.frame(maf.median=rep(0,length(tum)),
                      maf.IQR=rep(0,length(tum)),
                      type=rep("data", length(tum)))
  for(k in 1:length(tum)){
    ind<-which(dat.temp$tumor==tum[k])
    en<-dat.temp[ind,]
    data_M6[k,1:2]<-c(mean(en$median),
                      mean(en$IQR))
  }
  rm(dat.temp)
  
  return(list(data_M1_M2, data_M3, data_M4, data_M5, data_M6))
  
}

# ---------------------------------------------------------
# ---------------------------------------------------------
# ABC PRIORS
# ---------------------------------------------------------
# ---------------------------------------------------------

N.prior <- function(N_bounds){
  
  
  
  
  
  
}
