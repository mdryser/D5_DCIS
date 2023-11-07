#-----------------------------------------------
#-----------------------------------------------
#-- MAF-distributions across all DCIS spots
#-- Use for paper Figure 3B
#-- Updated: 09/13/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------


#-----------------------------------------------
#-- Load packages
#-----------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggExtra)

#-----------------------------------------------
#-- CHOICE
#-----------------------------------------------
## Parameter is TRUE if generating data_M6.csv: cap it at .5 (mutation there vs not)
## Parameter is FALSE if generating the paper figure (don't cap it)
modeling.version=FALSE


#-----------------------------------------------
#-- Preparation
#-----------------------------------------------

#-- Mutation present threshold
post.mut.thr<-.95
#-- Fraction of spots covered to assess public mutation
public.spot<-.9
#-- Only include spots with at least this many valid mutation reads
mut.min<-2


general_path <- "/Users/mdr30/Library/CloudStorage/Box-Box/0_D5_Paper/4_Post_processing/"
general_name_Mutations     <-  "Data/Cleaned_Mutation_calls/"
general_name_MAF   <- "Data/Cleaned_VAF/"


#-- Import the relevant files
Mut_files   <- list.files(path=general_name_Mutations, full.names = FALSE)
MAF_files  <- list.files(path=general_name_MAF, full.names = FALSE)


tumors<- sub(".*_", "", Mut_files)
tumors<- sub("csv.*", "", tumors)
tumors<-paste("Tumor_", tumors, sep="")

# Compiling all the MAFs
master.maf<-data.frame(tumor="a",
                       spot="b",
                       mean=0,
                       median=0,
                       sd=0,
                       IQR=0)

#-----------------------------------------------
#-- Main loop across tumors
#-----------------------------------------------

for(tum in tumors ){
  
  #-- Fetch the mutation data
  ind<-grep(tum ,Mut_files, fixed=TRUE, value=TRUE)
  dat.mut<-read.csv(file=paste(general_name_Mutations,ind,sep="/")) %>%
    filter(Spot!="Control") %>%
    filter(Histology=="DCIS") %>%
    remove_rownames() %>% 
    column_to_rownames(var = "Spot") %>%
    select(contains("chr"))
  
  #-- Fetch MAF files
  ind<-grep(tum ,MAF_files, fixed=TRUE, value=TRUE)
  dat.MAF<-read.csv(file=paste(general_name_MAF, ind,sep="/")) %>%
    filter(Spot!="Control") %>%
    filter(Histology=="DCIS") %>%
    remove_rownames() %>% 
    column_to_rownames(var = "Spot") %>%
    select(contains("chr"))
  
  #-- Only keep mutations that are present in at least one spot
  ind<-which(as.vector(apply(dat.mut, 2, function(x) any(x>post.mut.thr))))
  dat.mut<-dat.mut[,ind] 
  dat.MAF<-dat.MAF[,ind]
  rm(ind)
  
  #-- Check that the ordering is the same
  if(any(colnames(dat.mut)!=colnames(dat.MAF))){ warning('Issue with column or row ordering')}
  if(any(rownames(dat.mut)!=rownames(dat.MAF))){ warning('Issue with column or row ordering')}
  
  #-- Create a matrix MAF.redox where we determine the MAF of all mutations present at > post.mut.thr
  MAF.redox<-as.matrix(dat.MAF)
  ind<-which(dat.mut<post.mut.thr)
  MAF.redox[ind]<- NA
  
  #-- Only keep spots with at least mut.min mutations present
  valid.spots<-which(as.vector(rowSums(!!MAF.redox, na.rm = TRUE))>=mut.min) # number of valid reads in each spot
  
  MAF.redox<-MAF.redox[valid.spots,]
  
  ## HERE WATCH OUT
  ## ONLY for the modeling
  if(modeling.version==TRUE){
    MAF.redox<-pmin(MAF.redox, .5)
  }
  
  
  local.maf<-data.frame(tumor=rep(tum, length(valid.spots)),
                        spot=rownames(MAF.redox),
                        mean=as.vector(apply(MAF.redox, 1, function(x) mean(x, na.rm=TRUE))),
                        median=as.vector(apply(MAF.redox, 1, function(x) median(x, na.rm=TRUE))),
                        sd=as.vector(apply(MAF.redox, 1, function(x) sd(x, na.rm=TRUE))),
                        IQR=as.vector(apply(MAF.redox, 1, function(x) IQR(x, na.rm=TRUE))))
  
  #-- Concatenate
  master.maf<-rbind(master.maf,local.maf)
  
}

master.maf<-master.maf[-1,]


#-- Change the tumor name format to numeric label
master.maf$tumor<-sub("[.].*", "", master.maf$tumor)
master.maf$tumor<-as.numeric(as.character(sub(".*_", "", master.maf$tumor)))



#-----------------------------------------------
#-- OUTPUT SAVE & PLOTTING
#-----------------------------------------------

if(modeling.version){
  write.csv(master.maf, file="data_M6.csv", row.names = FALSE)
} else {
  
  write.csv(master.maf, file="Spot_MAF_summary_DCIS.csv", row.names = FALSE)
  
  
  p<-ggplot(master.maf, aes(x=median, y=IQR))+
    stat_density_2d(aes(fill = ..level..), geom = "polygon") +
    scale_fill_continuous(name="Density", low="lavenderblush", high="#CD534CFF") +
    geom_point(alpha=0.5, size = 1.5)+
    ylim(0,1)+xlim(0,1) +coord_fixed()+
    xlab("Variant allele frequency: median")+
    ylab("Variant allele frequency: IQR")+
    theme_bw()+
    theme(legend.position = c(1.12, .8))+
    theme(legend.text=element_text(size=10))+
    theme(legend.title = element_text(colour="black", size=10))+
    theme(text = element_text(size = 20)) 
  
  
  p<-ggExtra::ggMarginal(p, type = "density")
  
  
  #-- Save Figure 3B to PDF
  ggsave(filename = "Figure_3B.pdf",
         plot=p,
         width=20,
         height=20,
         units="cm")  
}

