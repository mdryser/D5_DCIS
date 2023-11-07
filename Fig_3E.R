#-----------------------------------------------
#-----------------------------------------------
#-- Spatial-genetic correlations: summary plot
#-- Use for paper Figure 3E
#-- Updated: 09/10/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------


#-----------------------------------------------
#-- Load packages
#-----------------------------------------------
library(readr)
library(ggplot2)
library(ggsci)
library(scales)
library(dplyr)


#-----------------------------------------------
#-- Initialize 
#-----------------------------------------------

#-- Color palette
pal1 <- pal_jco("default")(10)

#-- Seed
set.seed(1538)

#-- Tumor names to input to function
#-- Excluding 222:  only 2 DCIS spots
Tumornames <- c("03","17","24","53","65","66","88","91","168","173","191","211","213","225","227","286","315")

#-----------------------------------------------
#-- Parameters
#-----------------------------------------------

purity=.8 # tumor purity
cutoff=.05 # detection threshold
c0=.5 # prior prob of a mutation being absent (delta-weight)
p.mut.no<-c0*(1+pbeta(cutoff*purity, 1, 1/purity))
prob.present <- 1-p.mut.no

#-- Number of Monte Carlo realizations
n.sheets <- 1000

#-----------------------------------------------
#-- Main Function for Pearson R, given tumor
#-- Input: frame with spatial distances and mutation calls
#-----------------------------------------------

Find_Pearson <- function(sheet1){
  
  Tumor <- sheet1
  
  l<-dim(Tumor)[1]
  
  #-- Creating a table which contains each possible pairing of spots, spatial distance, and genetic distance
  GVS <- data.frame(matrix(nrow = l*(l-1)/2, ncol = 6))
  colnames(GVS) <- c("Spot_1", "Spot_2","Spatdist_abs", "Spatdist_norm", "Gendist_abs","Gendist_norm")
  
  #-- Create all pairwise combinations of spots
  GVS[,1:2]<-t(combn(Tumor$Spot, 2))
  
  
  #-- Spatial distance between each spot
  for(i in 1:dim(GVS)[1]){
    row1 <- Tumor[which(Tumor$Spot == GVS[i,1]),2:4]
    row2 <- Tumor[which(Tumor$Spot == GVS[i,2]),2:4]
    GVS$Spatdist_abs[i] <- dist(rbind(row1,row2), method = "euclidian")
  }
  
  #-- Filling in normalized spatial distance
  GVS$Spatdist_norm <-GVS$Spatdist_abs/max(GVS$Spatdist_abs)
  
  #-- Genetic distance between each spot
  for(i in 1:dim(GVS)[1]){
    row1 <- c(as.numeric(Tumor[which(Tumor$Spot == GVS[i,1]),13:dim(Tumor)[2]]))
    row2 <- c(as.numeric(Tumor[which(Tumor$Spot == GVS[i,2]),13:dim(Tumor)[2]]))
    GVS$Gendist_abs[i] <- dist(rbind(row1,row2), method = "manhattan")
  }
  
  #-- Normalized genetic distance
  GVS$Gendist_norm <- GVS$Gendist_abs/(dim(Tumor)[2]-12)
  
  #-- Pearson's R between spatial and genetic Data
  Pearson <- cor(GVS$Spatdist_norm,GVS$Gendist_norm,method = "pearson")
  
  #-- Returning Pearson's R value
  return(Pearson)
  
}


#-----------------------------------------------
#-- Function for Monte Carlo Sampling
#-----------------------------------------------

GVS_Rvalues <- function(DCIS_name, n.sheets, prob.present){
  
  #-- Reading in Tumor
  DCIS_path <- paste("Data/Cleaned_Mutation_calls/CleanedBayesian_Tumor_", DCIS_name, ".csv", sep = "")
  Tumor <- as.data.frame(read_csv(DCIS_path))
  Tumor <- Tumor[,-1]
  
  #-- Removing non-DCIS spots
  Tumor <- Tumor[which(Tumor$Histology == "DCIS"),]
  
  #------------------------------------
  #-- Replace LQT with prior probability 
  #------------------------------------
  
  # NA's denote a prob.present probability of the mutation being present
  for(i in 1:dim(Tumor)[1]){
    for(j in 13:dim(Tumor)[2]){
      if(is.na(Tumor[i,j])){
        Tumor[i,j] <- prob.present
      }
    }
  }
  
  #------------------------------------
  #-- Sample from posteriors
  #------------------------------------
  
  # Creating a list of length n.sheets with data frames of the size of the given tumor
  Allsheets <- vector(mode = "list", length = n.sheets)
  
  for(i in 1:n.sheets){
    
    binarysheet <- Tumor
    
    l<-dim(binarysheet)[2]
    
    #-- convert to a vector
    v<-c(data.matrix(binarysheet[,13:l]))
    
    binarysheet[,13:l]<-matrix(rbinom(length(v), 1, v), dim(binarysheet)[1], l-12)
    
    Allsheets[[i]] <- binarysheet
    
  }
  
  
  #------------------------------------
  #-- Store all the Pearson's R values
  #------------------------------------
  AllPearsons <- rep(0, n.sheets)
  
  for(i in 1:n.sheets){
    sheet1 <- Allsheets[[i]]
    AllPearsons[i] <- Find_Pearson(sheet1)
  }
  
  #-- Returning the plot and Pearson's R values
  return(AllPearsons)
  
}


#----------------------------------------
#----------------------------------------
#-- Main execution
#----------------------------------------
#----------------------------------------

GVS_store <- vector(mode = "list", length = length(Tumornames))
names(GVS_store) <- Tumornames

for(nam in Tumornames){
  
  print(nam)
  GVS_store[[nam]] <- GVS_Rvalues(nam,n.sheets,prob.present)
}

#-----------------------------
#-- Save the output
#-----------------------------

save(GVS_store, Tumornames, file="GVS_store.RData" )


#-----------------------------
#-- Load the output
#-----------------------------
load(file="GVS_store.RData")


#-----------------------------
#-- Compile Summary Statistics
#-----------------------------


#-- Building data table
Plot_data <- data.frame(matrix(nrow = length(Tumornames), ncol = 4))
colnames(Plot_data) <- c("Tumor","Mean","CI_low","CI_high")
Plot_data[,1] <- Tumornames

#-- Taking means
for(nam in Tumornames){
  Plot_data[which(Plot_data$Tumor == nam),"Mean"] <- mean(GVS_store[[nam]])
}

#-- Confidence interval 
for(nam in Tumornames){
  q<-quantile(GVS_store[[nam]], probs=c(.025, .975))
  
  Plot_data[which(Plot_data$Tumor == nam),"CI_low"] <- q[1]
  Plot_data[which(Plot_data$Tumor == nam),"CI_high"] <- q[2]
}

#-----------------------------
#-- ggplotting
#-----------------------------

#-- Insert tumor type
Plot_data$Status <- rep("Mixed")
for(i in 1:dim(Plot_data)[1]){
  if(Plot_data[i,"Tumor"] %in% c("03","17","24","211","213","222","225","227","315")){
    Plot_data[i,"Status"] <- "Pure"
  } 
}
Plot_data <- Plot_data[order(Plot_data$Status, decreasing = TRUE),]
Plot_data$Tumor <- factor(Plot_data$Tumor,levels = unique(Plot_data$Tumor))


ggplot(Plot_data,aes(x = Tumor, y = Mean, color = Status)) + 
  geom_point(size = 4) + 
  scale_color_manual(values= c("Pure" = pal1[3], "Mixed" = pal1[4]), labels = c("Mixed"="Sychronous DCIS & IBC","Pure"="DCIS")) +
  labs(y = "Pearson's R", x = "Tumor ID", title = "", color = "Tumor Status") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_errorbar(aes(ymin= CI_low, ymax= CI_high), width=.2, position=position_dodge(.9),size=1 ) +
  ylim(-1,1) +
  theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), plot.title = element_text(size = 18))+
  theme(axis.text.x = element_text(size=18),axis.text.y = element_text( size=18)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = pal1[[3]]) +
  theme(legend.position = "none")

#-----------------------------
#-- Save Figure 3E to PDF
#-----------------------------
ggsave("Figure_3E.pdf",
       width=25,
       height=10,
       units="cm")

#-----------------------------
#-- Write modeling output
#-----------------------------
write.csv(Plot_data, file="data_M3.csv")


#-----------------------------
#-- Statistical tests
#-----------------------------


#-- Test between pure and synchronous tumors
puretumors <- filter(Plot_data, Status == "Pure")
mixedtumors <- filter(Plot_data, Status == "Mixed")
t.test(puretumors$Mean, mixedtumors$Mean)
wilcox.test(puretumors$Mean, mixedtumors$Mean)


