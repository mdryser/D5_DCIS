#-----------------------------------------------
#-----------------------------------------------
#-- Expansion Index Summary
#-- Use for paper Figure 4D
#-- Updated: 09/13/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------


#-----------------------------------------------
#-- Load packages
#-----------------------------------------------

library(readr)
library(ggsci)
library(scales)
library(pracma)
library(dplyr)
library(ggplot2)

#-----------------------------------------------
#-- Load packages
#-----------------------------------------------

#-- Color palette
pal1 <- pal_jco("default")(10)

#-- Seed
set.seed(58901)

# Read in the tumors
# Exclude 222: only 2 DCIS spots
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


#-------------------------------------------------------
#-- Auxiliary function:  EI function value at x
#-------------------------------------------------------
EI.fun <- function(x, coord){

  ind<-which(coord$spots <= x)
  
  if(length(ind)){
    
    out<-max(coord$diam[ind])
    
  } else {
    
    out<-0
  }
  
  return(out)
}

#-- Vectorize the function 
EI.fun<-Vectorize(EI.fun, vectorize.args = "x")

#-------------------------------------------------------
#-- Function: EI for a given tumor input
#-------------------------------------------------------

EI_calc <- function(sheet){
  
  #-- Total number of DCIS spots
  spot.max<-dim(sheet)[1]
  
  #-- Maximum distance between any two DCIS spots
  diam.max<-max(dist(sheet[,2:4]))
  
  #-- Total number of mutations to check
  l<-dim(sheet)[2]-12
  
  #-- Keep record of spots and diameter for each mutation
  coord<-data.frame(spots=numeric(l),
                    diam=numeric(l))
  
  #-- loop through the mutations
  for(k in 1:l){
    
    #-- Find spots where mutation present
    ind<-which(sheet[,k+12]==1)
    
    #-- If there are at least two spots
    if(length(ind)>1){
      coord$diam[k]<-max(dist(sheet[ind,2:4]))/diam.max
      coord$spots[k]<-length(ind)/spot.max
    } else {   # Else fill with NA
      coord$diam[k]<-NA
      coord$spots[k]<-NA
    }
    
  }
  
  #-- Remove those mutations that were present in one or less spots
  coord<- na.omit(coord)
  EI_val<- integrate(EI.fun, 0, 1, coord=coord)
  return(EI_val$value)
  
}


#-------------------------------------------------------
#-- Monte Carlo sampling of Expansion Index calculation
#-------------------------------------------------------

EI_MC <- function(DCIS_name,n.sheets,prob.present){
  
  #-- Read in the tumor
  DCIS_path <- paste("Data/Cleaned_Mutation_calls/CleanedBayesian_Tumor_", DCIS_name, ".csv", sep = "")
  Tumor <- as.data.frame(read_csv(DCIS_path))
  Tumor <- Tumor[,-1]
  
  #-- Remove non-DCIS spots
  Tumor <- Tumor[which(Tumor$Histology == "DCIS"),]
  
  #-- Replacing Low-quality targets with prior probability

  for(i in 1:dim(Tumor)[1]){
    for(j in 13:dim(Tumor)[2]){
      if(is.na(Tumor[i,j])){
        Tumor[i,j] <- prob.present
      }
    }
  }
  
  #-- Sample from the posterior probability

  #-- Create list (length: n.sheets) of data frames of the size of the given tumor
  Allsheets <- vector(mode = "list", length = n.sheets)
  
  for(i in 1:n.sheets){
    binarysheet <- Tumor
    l<-dim(binarysheet)[2]
    v<-c(data.matrix(binarysheet[,13:l]))
    binarysheet[,13:l]<-matrix(rbinom(length(v), 1, v), dim(binarysheet)[1], l-12)
    
    Allsheets[[i]] <- binarysheet
    
  }
  
  # Extract the spot-diameter pairs (normalized) for each sheet
  EIvalue_list <- rep(NA, n.sheets)
  
  for(i in 1:n.sheets){
    EIvalue_list[i] <- EI_calc(Allsheets[[i]])
  }
  
  return(EIvalue_list)
  
}


#--------------------------------------
#-- MAIN LOOP
#--------------------------------------

# create the storage list
EI_store <- vector(mode = "list", length = length(Tumornames))
names(EI_store) <- Tumornames

# loop over the vectors
for(nam in Tumornames){
  print(nam)
  EI_store[[nam]] <- EI_MC(nam,n.sheets,prob.present)
}

#-- Save MC ouptut
save(EI_store, file="EI_store.RData")

#-- Reload MC output
load(file="EI_store.RData")

#--------------------------------------
#-- Summary Statistics
#--------------------------------------

#-- Building data table
Plot_data <- data.frame(matrix(nrow = length(Tumornames), ncol = 4))
colnames(Plot_data) <- c("Tumor","Mean","CI_low","CI_high")
Plot_data[,1] <- Tumornames

#-- Posterior means
for(nam in Tumornames){
  Plot_data[which(Plot_data$Tumor == nam),"Mean"] <- mean(EI_store[[nam]])
}

#-- Posterior prediction intervals
for(nam in Tumornames){
  q<-quantile(EI_store[[nam]], probs=c(.025, .975))
  
  Plot_data[which(Plot_data$Tumor == nam),"CI_low"] <- q[1]
  Plot_data[which(Plot_data$Tumor == nam),"CI_high"] <- q[2]
}

#-- Append tumor type
Plot_data$Status <- rep("Mixed")
for(i in 1:dim(Plot_data)[1]){
  if(Plot_data[i,"Tumor"] %in% c("03","17","24","211","213","222","225","227","315")){
    Plot_data[i,"Status"] <- "Pure"
  } 
}
Plot_data <- Plot_data[order(Plot_data$Status, decreasing = TRUE),]
Plot_data$Tumor <- factor(Plot_data$Tumor,levels = unique(Plot_data$Tumor))


#--------------------------------------
#-- Plotting
#--------------------------------------

ggplot(Plot_data,aes(x = Tumor, y = Mean, color = Status)) + 
  geom_point(size = 4) + 
  scale_color_manual(values= c("Pure" = pal1[3], "Mixed" = pal1[4]), labels = c("Mixed"="Sychronous DCIS & IBC","Pure"="DCIS")) +
  labs(y = "Expansion Index", x = "Tumor ID", title = "", color = "Tumor Status") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_errorbar(aes(ymin= CI_low, ymax= CI_high), width=.2, position=position_dodge(.9), size=1) +
  theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), plot.title = element_text(size = 18))+
  theme(axis.text.x = element_text(size=18),axis.text.y = element_text( size=18)) + 
  geom_hline(yintercept = 0.50, linetype = "dashed", color = pal1[[3]]) +
  coord_cartesian(ylim = c(0,1.0)) +
  theme( axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), plot.title = element_text(size = 18))+
  theme(axis.text.x = element_text(size=18),axis.text.y = element_text( size=18)) +
  theme(legend.position = "none")


ggsave("Figure_4D.pdf",
       width=25,
       height=10,
       units="cm")


#--------------------------------------
#-- Modeling output
#--------------------------------------

write.csv(Plot_data, file="data_M4.csv")


#--------------------------------
#-- Tumor-level statistics
#--------------------------------

#-- Pure vs synchronous tumors
puretumors <- filter(Plot_data, Status == "Pure")
mixedtumors <- filter(Plot_data, Status == "Mixed")

median(puretumors$Mean)
median(mixedtumors$Mean)
wilcox.test(puretumors$Mean, mixedtumors$Mean)
t.test(puretumors$Mean, mixedtumors$Mean)


