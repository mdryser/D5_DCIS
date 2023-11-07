#-----------------------------------------------
#-----------------------------------------------
#-- Tumor-level summary of final study cohort
#-- Use for Supplementary Table S1
#-- Updated: 09/13/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------

#-- Libraries
library(readr)
library(stringr)
library(dplyr)


#-- Selecting tumors
Tumornames <- c("03","17","24","211","213","222","225","227","315", "53","65","66","88","91","168","173","191","286")


#-- Creating an empty dataframe where all measures will be stored
Table1 <- data.frame(matrix(ncol = 12, nrow = length(Tumornames)))
row.names(Table1) <- Tumornames
colnames(Table1) <- c("Tumor","Slides","Total_Spots","Invasive_Spots","DCIS_Spots","Atypical_Spots","Non_Atypical_Spots","Normal_Spots", "Total_Mutations", "Total_reads","Low_Quality_Reads","Likely_Drivers")


#-----------------------------------------------
#- Main function that processes each tumor
#-----------------------------------------------

Tumor_data <- function(DCIS_name, Tab1){
  
  #-- importing data from the tumor 
  DCIS_path <- paste("Data/Cleaned_Mutation_calls/CleanedBayesian_Tumor_", DCIS_name, ".csv", sep = "")
  Tumor <- as.data.frame(read_csv(DCIS_path))
  Tumor <- Tumor[,-1]
  
  #-- Tumor name
  row <- which(rownames(Tab1) == DCIS_name)
  Tab1[row,1] <- DCIS_name
  
  #-- Slides and ordering
  slides <- c(unique(Tumor$Slide))
  slides <- slides[!is.na(slides)]
  order <- order(parse_number(slides))
  slides <- slides[c(order)]
  
  #-- Separating slide names in a single character string
  slides2 <- ""
  for(i in 1:length(slides)){
    slides2 <- paste(slides2,slides[i], sep = ",")
  } 
  slides2 <- substring(slides2,2)
  
  #-- Inserting slide names into table
  Tab1[row,2] <- slides2
  
  #-- Counting number of spots
  totalspots <- length.POSIXlt(Tumor) - 1
  Tab1[row,3] <- totalspots
  
  #-- Counting invasive spots
  invasivespots <- length(which(Tumor$Histology == "Invasive"))
  Tab1[row,4] <- invasivespots
  
  #-- Counting DCIS spots
  Dcisspots <- length(which(Tumor$Histology == "DCIS"))
  Tab1[row,5] <- Dcisspots
  
  #-- Non atypical spots
  Nonatypicalspots <- length(which(Tumor$Histology == "Non-Atypical"))  
  Tab1[row,6] <- Nonatypicalspots
  
  #-- Atypical spots
  Atypicalspots <- length(which(Tumor$Histology == "Atypical"))
  Tab1[row,7] <- Atypicalspots
  
  #-- Normal spots
  normalspots <- length(which(Tumor$Histology == "Normal"))
  Tab1[row,8] <- normalspots
  
  #-- Total Mutations
  totalmuts <- length(Tumor) - 12
  Tab1[row,9] <- totalmuts
  
  #-- Total reads
  Tab1[row,10] <- (totalmuts * totalspots)
  
  #-- Low quality reads
  numlqr <- 0
  for(i in 2:length.POSIXlt(Tumor)){
    for(j in 13:length(Tumor)){
      if(is.na(Tumor[i,j])){
        numlqr <- numlqr +1
      }
    }
  }
  Tab1[row,11] <- numlqr
  
  
  #-- Importing  driver status
  Mutations <- as.data.frame(read_csv("Data/Annotations_SNV/SNV_annotations_complete.csv"))  %>%
    mutate(Driver=replace(Driver, Driver=="", "intergenic")) %>%
    rename(Tumor=DCISsample)
  
  #-- Renaming NA chromosomes as X
  for(i in 1:length.POSIXlt(Mutations)){
    if(is.na(Mutations[i,"CHROM"])){
      Mutations[i,"CHROM"] <- "X"
    }
  }
  
  #-- Naming mutations and taking necessary rows
  Mutations$MutName <- c(rep(0,times = length.POSIXlt(Mutations)))
  for(i in 1:length.POSIXlt(Mutations)){
    Mutations[i,"MutName"] <- paste("chr",Mutations[i,"CHROM"],":",Mutations[i,"POS"],sep = "")
  }

  Mutations <- Mutations[,c("MutName","Tumor","Driver")]
  
  #-- Renaming Tumor names and NAs
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
  
  
  correctmuts <- Mutations[which(Mutations$Tumor == DCIS_name),]
  numdriver <- length(which(correctmuts[,"Driver"] == "likely"))
  
  #-- Completing table
  Tab1[row,12] <- numdriver
  
  return(Tab1)
}


#---------------------------------
#- Loop over all tumors
#---------------------------------

for(tum in Tumornames){
  Table1 <- Tumor_data(tum, Table1)
}


#---------------------------------
#- Final curation and column sums
#---------------------------------

Table1 <- Table1 %>%
  mutate(Slide.n=lengths(gregexpr(",", Slides)) + 1) %>%
  mutate(Benign_Spots=Non_Atypical_Spots+Atypical_Spots) %>%
  select(Slide.n, Total_Spots, DCIS_Spots, Invasive_Spots, Benign_Spots, Normal_Spots, Total_Mutations, Likely_Drivers, Total_reads, Low_Quality_Reads)

Table1 <- rbind(Table1, colSums(Table1))

#---------------------------------
#- Printout for Table S1
#---------------------------------

path <- "/Table_S1.csv"
write.csv(Table1, path)



