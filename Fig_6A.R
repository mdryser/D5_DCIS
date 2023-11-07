#-----------------------------------------------
#-----------------------------------------------
#-- Mutation Network Graphs
#-- Use for paper Figure 6A, 6D
#-- Updated: 09/13/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------

#-----------------------------------------------
#-- Load packages
#-----------------------------------------------
library(readr)
library(plotly)
library(ggsci)
library(scales)
library(stringr)

#-----------------------------------------------
#-- Preparations
#-----------------------------------------------

#-- Color palette
pal1 <- pal_jco("default")(10)

#-- Read in the tumor of interest
DCIS_name <- "222"
DCIS_path <- paste("Data/Cleaned_Mutation_calls/CleanedBayesian_Tumor_", DCIS_name, ".csv", sep = "")
Tumor <- as.data.frame(read_csv(DCIS_path))
Tumor$X1 <- NULL
colnames(Tumor) <-  gsub(" ", "", colnames(Tumor), fixed = TRUE)
Tumor <- Tumor[-1,]
row.names(Tumor) <- 1:length.POSIXlt(Tumor)



#-----------------------------------------------
#-- Mutation calls
#-----------------------------------------------
for(row in 1:length.POSIXlt(Tumor)){
  for(col in 14:length(Tumor)){
    if(is.na(Tumor[row,col])){
      Tumor[row,col] <- 0
    }
    else if(Tumor[row,col] >= 0.95){
      Tumor[row,col] <- 1
    }
    else if(Tumor[row,col] <= 0.05){
      Tumor[row,col] <- 0
    } 
    else(Tumor[row,col] <- 0)
  }
}

#-----------------------------------------------
#-- Spot histologies: 4 categories
#-----------------------------------------------
for(i in 1:length.POSIXlt(Tumor)){
  if(Tumor[i,"Histology"] == "Atypical" | Tumor[i,"Histology"] == "Non-Atypical"){
    Tumor[i,"Histology"] <- "Benign"
  }
}


#-----------------------------------------------
#-- Data tables
#-----------------------------------------------
Nodes <- data.frame(matrix(ncol = 6, nrow = length.POSIXlt(Tumor)))
colnames(Nodes) <- c("Spot_name","group","color","x","y","ID")

Links <- data.frame(matrix(ncol = 7, nrow = length.POSIXlt(Tumor)**2))
colnames(Links) <- c("Source","Target","Weight","Group", "IDsource","IDtarget","Color")

#-----------------------------------------------
#-- Fill in nodes
#-----------------------------------------------
Nodes$Spot_name <- Tumor$Spot
Nodes$group <- Tumor$Histology
Nodes$color <- Tumor$Slide

#-- Filling in x position
for(i in 1:length.POSIXlt(Nodes)){
  if(Nodes[i,"group"] == "Normal"){
    Nodes[i,"x"] <- 0.05
  }else if(Nodes[i,"group"] == "Benign"){
    Nodes[i,"x"] <- 0.33
  } else if(Nodes[i,"group"] == "DCIS"){
    Nodes[i,"x"] <- 0.66
  } else if(Nodes[i,"group"] == "Invasive"){
    Nodes[i,"x"] <- 0.95
  }
}

#-- Filling in y position based on slide
slide_names <- data.frame(unique(Tumor$Slide))
slide_names$numbers <- parse_number(slide_names$unique.Tumor.Slide.)
slide_names <- slide_names[order(slide_names$numbers),]
slide_names <- c(slide_names$unique.Tumor.Slide.)

slideinterval <- 1/length(slide_names)
for(i in 1:length.POSIXlt(Nodes)){
  if(Nodes[i,"color"] == slide_names[1]){
    Nodes[i,"y"] <- runif(1, min = 0, max = 0 + slideinterval)
  } else if(Nodes[i,"color"] == slide_names[2]){
    Nodes[i,"y"] <- runif(1, min = slideinterval, max = 2*slideinterval)
  } else if(length(slide_names) >= 3 & Nodes[i,"color"] == slide_names[3]){
    Nodes[i,"y"] <- runif(1, min = 2*slideinterval, max = 3*slideinterval)
  } else if(length(slide_names) >= 4 & Nodes[i,"color"] == slide_names[4]){
    Nodes[i,"y"] <- runif(1, min = 3*slideinterval, max = 4*slideinterval)
  } else if(length(slide_names) >= 5 & Nodes[i,"color"] == slide_names[5]){
    Nodes[i,"y"] <- runif(1, min = 4*slideinterval, max = 5*slideinterval)
  }
}


#-- Filling in color based on slide
for(i in 1:length.POSIXlt(Nodes)){
  if(Nodes[i,"color"] == slide_names[1]){
    Nodes[i,"color"] <- pal1[1]
  } else if(Nodes[i,"color"] == slide_names[2]){
    Nodes[i,"color"] <- pal1[2]
  } else if(length(slide_names) >= 3 & Nodes[i,"color"] == slide_names[3]){
    Nodes[i,"color"] <- pal1[7]
  } else if(length(slide_names) >= 4 & Nodes[i,"color"] == slide_names[4]){
    Nodes[i,"color"] <- pal1[8]
  } else if(length(slide_names) >= 5 & Nodes[i,"color"] == slide_names[5]){
    Nodes[i,"color"] <- pal1[6]
  }
}

#- Filling in ID
Nodes$ID <- row.names(Nodes)



#-----------------------------------------------
#-- Driver/passenger status
#-----------------------------------------------
#-- Import mutation driver data and extracting necessary data
Mutations <- as.data.frame(read_csv("Data/Annotations_SNV/SNV_annotations_complete.csv"))  %>%
  mutate(Driver=replace(Driver, Driver=="", "intergenic")) %>%
  rename(Tumor=DCISsample)

#-- Rename NA chromosomes as X
for(i in 1:length.POSIXlt(Mutations)){
  if(is.na(Mutations[i,"CHROM"])){
    Mutations[i,"CHROM"] <- "X"
  }
}

#-- Name mutations and taking necessary rows
Mutations$MutName <- c(rep(0,times = length.POSIXlt(Mutations)))
for(i in 1:length.POSIXlt(Mutations)){
  Mutations[i,"MutName"] <- paste("chr",Mutations[i,"CHROM"],":",Mutations[i,"POS"],sep = "")
}

Mutations <- Mutations[,c("MutName","Tumor","Driver")]

#-- Rename Tumor names and NAs
for(i in 1:length.POSIXlt(Mutations)){
  Mutations[i,"Tumor"] <- as.character(parse_number(Mutations[i,"Tumor"]))
  if(is.na(Mutations[i,"Driver"])){
    Mutations[i,"Driver"] <- "unlikely"
  }
}

for(i in 1:length.POSIXlt(Mutations)){
  if(Mutations[i,"Tumor"] == "3"){
    Mutations[i,"Tumor"] <- "03"
  }
}

#-- Function finding status of given tumors mutations
Tumor_muts <- Mutations[which(Mutations[,"Tumor"] == DCIS_name),]


#-- Fill in the links
Links$Source <- c(rep(Nodes$Spot_name, each = length.POSIXlt(Nodes)))
Links$Target <- c(rep(Nodes$Spot_name, times = length.POSIXlt(Nodes)))

#-- Remove rows with duplicated pairing or with 2 of the same spot
remove <- c()
for(i in 1:length.POSIXlt(Links)){
  if(Links[i,1] >= Links[i,2]){
    remove <- c(remove,i)
  }
}
if(length(remove)){Links <- Links[-remove,]}


#-- Find number of shared mutations for the weight
for(i in 1:length.POSIXlt(Links)){
  row1 <- c(as.numeric(Tumor[which(Tumor[,"Spot"] == Links[i,"Source"]),14:length(Tumor)]))
  row2 <- c(as.numeric(Tumor[which(Tumor[,"Spot"] == Links[i,"Target"]),14:length(Tumor)]))
  Links[i,"Weight"] <- length(which(row1[which(row1 == row2)] == 1))
}
Links$Weight <- (Links$Weight/(length(Tumor)-12)) * 100


#-- Find status of connections and removing those in the same layer
for(i in 1:length.POSIXlt(Links)){
  Hist1 <- Nodes[which(Nodes[,"Spot_name"] == Links[i,"Source"]),"group"]
  Hist2 <- Nodes[which(Nodes[,"Spot_name"] == Links[i,"Target"]),"group"]
  
#-- If same layer: remove
  if(Hist1 == Hist2){
    Links[i,"Group"] <- "Same"
  } else if(Hist1 == "Benign" & Hist2 == "Invasive"){
    Links[i,"Group"] <- "Skip"
  } else if(Hist1 == "Invasive" & Hist2 == "Benign"){
    Links[i,"Group"] <- "Skip"
  } else{
    Links[i,"Group"] <- "Different"
  }
}

Links <- Links[which(Links[,"Group"] == "Different"),]


#-- Rearrange order of connection
for(i in 1:length.POSIXlt(Links)){
  Hist1 <- Nodes[which(Nodes[,"Spot_name"] == Links[i,"Source"]),"group"]
  Hist2 <- Nodes[which(Nodes[,"Spot_name"] == Links[i,"Target"]),"group"]
  
  sourcenum <- Links[i,"Source"]
  targnum <- Links[i,"Target"]
  
  if(Hist1 == "Invasive" & Hist2 == "DCIS"){
    Links[i,"Source"] <- targnum
    Links[i,"Target"] <- sourcenum
  } else if(Hist1 == "Invasive" & Hist2 == "Benign"){
    Links[i,"Source"] <- targnum
    Links[i,"Target"] <- sourcenum
  } else if(Hist1 == "DCIS" & Hist2 == "Benign"){
    Links[i,"Source"] <- targnum
    Links[i,"Target"] <- sourcenum
  } else if(Hist1 == "Invasive" & Hist2 == "Normal"){
    Links[i,"Source"] <- targnum
    Links[i,"Target"] <- sourcenum
  } else if(Hist1 == "DCIS" & Hist2 == "Normal"){
    Links[i,"Source"] <- targnum
    Links[i,"Target"] <- sourcenum
  } else if(Hist1 == "Benign" & Hist2 == "Normal"){
    Links[i,"Source"] <- targnum
    Links[i,"Target"] <- sourcenum
  }
}

#-- Fill in the color based on driver status

#-- Make all connections grey
Links$Color <- pal1[3]

#-- Rename tumor columns
for(mutation in colnames(Tumor)[14:length(Tumor)]){
  mut_position <- str_split(mutation,":")
  mut_position <- paste(mut_position[[1]][1],":",mut_position[[1]][2], sep = "")
  colnames(Tumor)[which(colnames(Tumor) == mutation)] <- mut_position
}

for(mut in Tumor_muts$MutName){
  if(Tumor_muts[which(Tumor_muts[,"MutName"] == mut),"Driver"] == "likely"){
    for(i in 1:length.POSIXlt(Links)){
      sourcerow <- Tumor[which(Tumor[,"Spot"] == Links[i,"Source"]),mut]
      targetrow <- Tumor[which(Tumor[,"Spot"] == Links[i,"Target"]),mut]
      if(sourcerow == 1 & targetrow == 1){
        Links[i,"Color"] <- pal1[9]
      }
    }
  }
}


#-- Fill in ID source and target
for(i in 1:length.POSIXlt(Links)){
  Links[i,"IDsource"] <- Nodes[which(Nodes[,"Spot_name"] == Links[i,"Source"]),"ID"] 
  Links[i,"IDtarget"] <- Nodes[which(Nodes[,"Spot_name"] == Links[i,"Target"]),"ID"] 
}
for(i in 1:length.POSIXlt(Links)){
  Links[i,"IDsource"] <- as.numeric(Links[i,"IDsource"]) -1
  Links[i,"IDtarget"] <- as.numeric(Links[i,"IDtarget"]) -1
}


#-- Remove links from benign to invasive
remove <- c()
for(i in 1:length.POSIXlt(Links)){
  Hist1 <- Nodes[which(Nodes[,"Spot_name"] == Links[i,"Source"]),"group"]
  Hist2 <- Nodes[which(Nodes[,"Spot_name"] == Links[i,"Target"]),"group"]
  if(Hist1 == "Benign" & Hist2 == "Invasive"){
    remove <- c(remove,i)
  } else if(Hist1 == "Normal" & Hist2 == "Invasive"){
    remove <- c(remove,i)
  } else if(Hist1 == "Normal" & Hist2 == "DCIS"){
    remove <- c(remove,i)
  } 
}
if(length(remove)){Links <- Links[-remove,]}

#-- Links with 0 weight won't show up
mycol <- rgb(0, 0, 255, max = 255, alpha = 0, names = "blue50")
Links[which(Links$Weight == 0), "Color"] <- mycol

#-----------------------------------------------
#-- The Final Plots
#-----------------------------------------------
p <- plot_ly(
  type = "sankey",
  domain = c(
    x =  c(0,1),
    y =  c(0,1)
  ),
  valueformat = "",
  valuesuffix = "% Mutations Shared",
  
  #-- Add nodes
  node = list(
    label = Nodes$Spot_name,
    color = Nodes$color,
    x = Nodes$x,
    y = Nodes$y,
    thickness = 10,
    length = 10,
    line = list(
      color = "black",
      width = 0.1
    )
  ),
  
  #-- Add links
  link = list(
    source = Links$IDsource,
    target = Links$IDtarget,
    value =  rep(1, times = length.POSIXlt(Links)),
    #value =  Links$Weight,
    color =  Links$Color,
    label =  Links$Group
    
  )
)

p <- p %>% layout(
  title = paste ("Driver Mutation Flow Tumor ", DCIS_name),
  font = list(
    size = 10,
    color = 'black'
  ),
  xaxis = list(showgrid = F, zeroline = F, showticklabels = F),
  yaxis = list(showgrid = F, zeroline = F, showticklabels = F),
  plot_bgcolor = 'white',
  paper_bgcolor = 'white'
)

p






