#-----------------------------------------------
#-----------------------------------------------
#-- Spatially clustered DCIS mutation heatmaps
#-- Use for paper Figure 2C, 2D
#-- Updated: 09/12/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------


#-----------------------------------------------
#-- Load packages
#-----------------------------------------------
library(gplots)
library(dplyr)
library(Rtsne)
library(ggsci)
library(readr)
library(scales)
library(stringr)
library(ggpubr)


#-----------------------------------------------
#-- Initialize stuff
#-----------------------------------------------
set.seed(223)

#-- Setting working directory and palette
pal1 <- pal_jco("default")(10)


heatmap_table <- function(DCIS_name){
  
  #-- Reading in Tumor
  DCIS_path <- paste("Data/Cleaned_Mutation_calls/CleanedBayesian_Tumor_", DCIS_name, ".csv", sep = "")
  Tumor <- as.data.frame(read_csv(DCIS_path))
  Tumor <- Tumor[-1,]
  
  dat.mut <- Tumor[which(Tumor$Histology == "DCIS"),]
  #-- dat.mut <- Tumor
  
  # All NA's denote a 0.5 probability of the mutation being present
  for(i in 1:length.POSIXlt(dat.mut)){
    for(j in 13:length(dat.mut)){
      if(is.na(dat.mut[i,j])){
        dat.mut[i,j] <- 0.5
      }
    }
  }
  
  #-- Do the spatial alignment, except for DCIS-222
  
  if(DCIS_name!="222"){
    spat<- dat.mut %>% select(x.pos, y.pos, z.pos)
    tsne_out <- Rtsne(as.matrix(spat),pca=TRUE,dims = 1,perplexity=1,theta=0) # Run TSNE
    dat.mut <- dat.mut[order(tsne_out$Y),]
  }
  
  row.names(dat.mut) <- dat.mut$Spot
  dat.mut <- as.data.frame(dat.mut)
  
  #-- Manipulating data table for input to heatmaps
  dat.mut <- dat.mut %>%
    select(which(grepl("chr", colnames(dat.mut), fixed=TRUE))) %>%
    mutate_all(as.numeric)
  dat.mut<-t(as.matrix(dat.mut))
  
  #-- Turning values into categorical (0, 1, and .5=NA)
  dat.mut[dat.mut <.05]<-0
  dat.mut[dat.mut >.95]<-1
  dat.mut[dat.mut >0 & dat.mut<1]<-.5
  dat.mut <- as.data.frame(dat.mut)
  
  
  #-- Changing spot histology to 3 categories
  for(i in 1:length.POSIXlt(Tumor)){
    if(Tumor[i,"Histology"] == "Atypical" | Tumor[i,"Histology"] == "Non-Atypical" | Tumor[i,"Histology"] == "Normal"){
      Tumor[i,"Histology"] <- "Benign"
    }
  }
  
  #-- Making a vector to give spot colors
  spotcolors <- data.frame(matrix(ncol = 2, nrow = length(dat.mut)))
  colnames(spotcolors) <- c("Spot","Color")
  spotcolors$Spot <- colnames(dat.mut)
  for(i in 1:length.POSIXlt(spotcolors)){
    spotcolors[i,"Color"] <- Tumor[which(Tumor[,"Spot"] == spotcolors[i,"Spot"]),"Histology"]
  }
  for(i in 1:length.POSIXlt(spotcolors)){
    if(spotcolors[i,"Color"] == "Benign"){
      spotcolors[i,"Color"] <- pal1[7]
    } else if(spotcolors[i,"Color"] == "DCIS"){
      spotcolors[i,"Color"] <- pal1[5]
    } else if(spotcolors[i,"Color"] == "Invasive"){
      spotcolors[i,"Color"] <- pal1[4]
    }
  }
  spotcolors <- c(spotcolors$Color)
  
  #-- Changing mutation names
  namevector <- c()
  for(i in 1:length.POSIXlt(dat.mut)){
    indname <- str_split(row.names(dat.mut)[i],":")
    indname <- indname[[1]][3]
    namevector <- c(namevector,indname)
  }
  row.names(dat.mut) <- make.names(namevector, unique = TRUE)
  
  
  return(list(dat.mut,spotcolors))
  
}

#-- Making a new plot to get legend from
forplotting <- data.frame(matrix(ncol = 3,nrow = 3))
colnames(forplotting) <- c("x.pos","y.pos","Histology")
forplotting$x.pos <- c(0.5,1,1.5)
forplotting$y.pos <- c(0.5,1,1.5)
forplotting$Histology <- c("Benign","DCIS","Invasive")
exampleplot <- ggplot(forplotting,aes(x = x.pos, y = y.pos)) +
  geom_point(size = 5, aes(color = Histology)) +
  scale_colour_manual(values = c("Benign" = pal1[7], "DCIS" = pal1[5],"Invasive" = pal1[4]), name = "Spot Histology", drop = FALSE) +
  theme(legend.position = "top", legend.box = "vertical") +
  theme_classic()
plotlegend <- get_legend(exampleplot)



#-----------------------------------------------
#-- Loop over all tumors
#-----------------------------------------------
tum.nam.vec<-c("03", "17", "24", "211", "222", "213",  "225", "227", "315", "53", "65", "66", "88", "91", "168", "173", "191", "286")

for(Tumorname in tum.nam.vec){

  #-- Step 1: Call the pdf command to start the plot
  pdf(file = paste("Figure_", Tumorname, ".pdf", sep=""),   # The directory you want to save the file in
      width = 10, # The width of the plot in inches
      height = 6) # The height of the plot in inches
  
  output <- heatmap_table(Tumorname)
  X1 <- as.matrix(output[[1]])
  spotcolors <- output[[2]]
  heatmap.2(X1,
            distfun = function (y) dist(y,method = "manhattan"),
            hclustfun = hclust, 
            scale = "none",
            col = pal1[c(1,3,2)],
            key=FALSE,
            tracecol=NA,
            key.title="Mutation Key",
            # colCol = spotcolors,
            revC = TRUE,
            Rowv=TRUE,
            Colv=FALSE,
            trace="none",
            cexRow=.6, 
            main="",
            margins = c(4,8),
            density.info = "none",
            breaks = c(0,0.05,0.95,1),
            key.xlab = "Probability",
            dendrogram = "none")
  
  
  dev.off()
  
}


