#-----------------------------------------------
#-----------------------------------------------
#-- Multiclonal Invasion
#-- Use for paper Figure 6H
#-- Updated: 09/13/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------


#-----------------------------------------------
#-- Load packages
#-----------------------------------------------

#-----------------------------------------------
#-- Load packages
#-----------------------------------------------
library(dplyr)
library(ggplot2)
library(readr)
library(ggsci)


#-----------------------------------------------
#-- Preparations
#-----------------------------------------------

#-- Color palette
pal1 <- pal_jco("default")(10)

#-- Seed
set.seed(58901)

#-- Tumor names for this analysis
Tumornames <- c("53", "65", "66", "88", "91", "168", "173", "286")



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

#-- Percentage of spots a public mutation must occupy to be public
thresh.public <- 0.90

#-----------------------------------------------
#-- Multiclonal invasion function
#-----------------------------------------------
Multi_Invasion <- function(prob.present,n.sheets,thresh.public, DCIS_name){
  
  #-- Reading in tumor data
  DCIS_path <- paste("Data/Cleaned_Mutation_calls/CleanedBayesian_Tumor_", DCIS_name, ".csv", sep = "")
  Tumor <- as.data.frame(read_csv(DCIS_path))
  Tumor <- Tumor[-1,-1]
  
  #-- At least 2 DCIS and invasive spots each
  if(length(which(Tumor[,"Histology"] == "Invasive")) < 2){
    if(length(which(Tumor[,"Histology"] == "DCIS")) < 2){
      return(NA)
    }
  }
  
  #-- NA values replaced by prior probability
  for(i in 1:length.POSIXlt(Tumor)){
    for(j in 13:length(Tumor)){
      if(is.na(Tumor[i,j])){
        Tumor[i,j] <- prob.present
      }
    }
  }
  
  
  #-----------------------------------------------
  #-- Sampling from posterior probabilities
  #-----------------------------------------------
  
  #-- List of length n.sheets with dataframes of the size of the given tumor
  Allsheets <- vector(mode = "list", length = n.sheets)
  
  for(i in 1:n.sheets){
    binarysheet <- Tumor
    for(row in 1:dim(binarysheet)[1]){
      for(column in 13:dim(binarysheet)[2]){
        binarysheet[row,column] <- sample(c(0,1),1,prob = c(1-binarysheet[row,column],binarysheet[row,column]))
      }
    }
    
    Allsheets[[i]] <- binarysheet
  }
  
  #-- Table to keep track of # of private invading mutations
  private_invaders <- data.frame(matrix(ncol = 3, nrow = n.sheets))
  colnames(private_invaders) <- c("Total_Mutations","Private_Invading_Mutations","Distinct_Invasional_Patterns")
  private_invaders$Total_Mutations <- c(rep(dim(Tumor)[2]-12))
  
  #-- Check each mutation for invasion
  for(i in 1:length(Allsheets)){
    currentsheet <- Allsheets[[i]]
    possDCIS <- round(length(which(currentsheet[,"Histology"] == "DCIS")) * thresh.public)
    possInvasive <- round(length(which(currentsheet[,"Histology"] == "Invasive")) * thresh.public)
    invasion_counter <- 0
    whichmuts <- c()
    for(col in colnames(currentsheet)[13:length(currentsheet)]){
      presentspots <- currentsheet[which(currentsheet[,col] == 1),]
      if(length(which(presentspots[,"Histology"] == "DCIS")) >= 1 & length(which(presentspots[,"Histology"] == "DCIS")) < possDCIS){
        if(length(which(presentspots[,"Histology"] == "Invasive")) >= 1 & length(which(presentspots[,"Histology"] == "Invasive")) < possInvasive){
          invasion_counter <- invasion_counter + 1
          whichmuts <- c(whichmuts,col)
        }
      }
    }
    private_invaders[i,2] <- invasion_counter
    
    #-- Find distinct patterns
    currentsheet <- currentsheet[(currentsheet$Histology == "DCIS" | currentsheet$Histology == "Invasive"),whichmuts]
    currentsheet <- as.data.frame(t(currentsheet))
    private_invaders[i,3]<-dim(unique(currentsheet))[1]
    
  }
  
  return(private_invaders)
  
}


#-----------------------------------------------
#-- Loop over all tumors and do MC
#-----------------------------------------------
Tumor_Invaders <- vector(mode = "list", length = length(Tumornames))
names(Tumor_Invaders) <- Tumornames
for(nam in Tumornames){
  print(nam)
  Tumor_Invaders[[nam]] <- Multi_Invasion(prob.present, n.sheets, thresh.public, nam)
}

#-- Creating plot data
PlotData <- data.frame(matrix(ncol = 5, nrow = length(Tumor_Invaders)))
colnames(PlotData) <- c("Tumor","Average_Private_Invading_Mutations", "Average_Unique_PIM", "PI_Low","PI_High")
PlotData$Tumor <- names(Tumor_Invaders)

for(i in 1:dim(PlotData)[1]){
  getdata <- Tumor_Invaders[[PlotData[i,"Tumor"]]]
  
  #-- Mean
  PlotData$Average_Private_Invading_Mutations[i] <- mean(getdata$Private_Invading_Mutations)
  PlotData$Average_Unique_PIM[i] <- mean(getdata$Distinct_Invasional_Patterns)
  
  #-- 95% interval
  huh<-quantile(getdata$Distinct_Invasional_Patterns, probs = c(.025, .975))
  PlotData$PI_Low[i] <- huh[1]
  PlotData$PI_High[i] <- huh[2]
}

#-----------------------------------------------
#-- Save the summary data
#-----------------------------------------------
save(PlotData, file="multiclonal_invasion_1000.RData")

#-----------------------------------------------
#-- Load the summary data
#-----------------------------------------------
load(file="multiclonal_invasion_1000.RData")


#-----------------------------------------------
#-- Plot the summary
#-----------------------------------------------
PlotData$Tumor <- factor(PlotData$Tumor, levels = PlotData$Tumor)
ggplot(PlotData, aes(x = Tumor, y = Average_Unique_PIM)) + 
  geom_col(fill = pal1[3], alpha=.6) +
  labs(y = "Multiclonal invasion, n", title = "", x="Tumor ID") +
  geom_errorbar(aes(ymin= PI_Low, ymax= PI_High), width=.2, position=position_dodge(.9)) + 
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(text = element_text(size = 20))  +
  theme( # remove the vertical grid lines
    panel.grid.major.x = element_blank() ,
    # explicitly set the horizontal lines (or they will disappear too)
    panel.grid.major.y = element_blank() 
  )

#-- Save Figure 6H to PDF
ggsave(file="Figure_6H.pdf",
       height=12,
       width=20,
       units="cm")

