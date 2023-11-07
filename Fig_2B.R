#-----------------------------------------------
#-----------------------------------------------
#-- Circos plot for D5 summary
#-- Use for paper Figure 2B
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
library(circlize) 
library(gplots)
library(dplyr)
library(Rtsne)
library(ggsci)
library(readr)
library(scales)
library(stringr)
library(ggpubr)
library(ComplexHeatmap)

#-----------------------------
#-- Initialize 
#-----------------------------

#-- Color palette
pal1 <- pal_jco("default")(10)

#-- Mutation present threshold
post.mut.thr<-.95

#-- Fraction of spots covered to assess public mutation
public.spot<-.9

#-- Path name
general_name_Mutations <- "Data/Cleaned_Mutation_calls/"

#-- Import driver mutations
driver.tab<-read.csv(file="Data/Annotations_SNV/Driver_mutations_All.csv")

#-- Import the mutation files
Mut_files   <- list.files(path=general_name_Mutations, full.names = FALSE)

#-- Process the mutation files
tumors<- sub(".*_", "", Mut_files)
tumors<- sub("csv.*", "", tumors)
tumors<-paste("Tumor_", tumors, sep="")

pure.tum<-c("Tumor_03.",
            "Tumor_17.",
            "Tumor_24.",
            "Tumor_211.",
            "Tumor_213.",
            "Tumor_222." ,
            "Tumor_225.",
            "Tumor_227.",
            "Tumor_315." )

mixed.tum<- c("Tumor_53.",
              "Tumor_65.", 
              "Tumor_66.", 
              "Tumor_88.", 
              "Tumor_91.",
              "Tumor_168.",
              "Tumor_173.",
              "Tumor_191.",
              "Tumor_286.")

tumors<-c(pure.tum, mixed.tum)


#-- Initialize the final matrix

df <- data.frame(matrix(ncol = 105, nrow = 0))
colnames(df)[1:5] <- c("Tumor","Spot", "Type", "Histology", "Driver")
colnames(df)[6:105]<- paste("mut", c(1:100))


#-----------------------------
#-- Main Loop across Tumors
#-----------------------------

for(tum in tumors ){
  
  #-- Fetch the mutation data
  ind<-grep(tum ,Mut_files, fixed=TRUE, value=TRUE)
  
  dat.mut<-read.csv(file=paste(general_name_Mutations,ind,sep="/")) %>%
    filter(Spot!="Control") %>%
    remove_rownames() %>% 
    select(-X)
  
  #-- Local dataframe
  df_loc <- data.frame(matrix(ncol = 105, nrow = dim(dat.mut)[1]))
  colnames(df_loc)[1:5] <- c("Tumor","Spot", "Type", "Histology", "Driver")
  colnames(df_loc)[6:105]<- paste("mut", c(1:100))
  
  
  #-- Fill in the first few columns
  df_loc$Tumor<-tum 
  df_loc$Spot<-dat.mut$Spot
  df_loc$Type<- ifelse(tum %in% pure.tum, "Pure", "Mixed" )
  df_loc$Histology <- dat.mut$Histology
  
  #-- Extract the mutation portion only
  dat.mut<-read.csv(file=paste(general_name_Mutations,ind,sep="/")) %>%
    filter(Spot!="Control") %>%
    remove_rownames() %>% 
    column_to_rownames(var = "Spot") %>%
    select(contains("chr"))
  
  dat.mut[dat.mut>post.mut.thr]<-1
  dat.mut[dat.mut<(1-post.mut.thr)]<-0
  dat.mut[dat.mut>(1-post.mut.thr) & dat.mut<post.mut.thr]<-.5
  dat.mut[is.na(dat.mut)]<-.5
  
  #-- Cluster by columns
  column_od<-hclust(dist(t(dat.mut),method = "manhattan"))$order
  dat.mut<-dat.mut[,column_od]
  
  #-- Identify the driver mutations
  mut.vec<-colnames(dat.mut)
  mut.vec.driver<-rep(0, length(mut.vec))
  
  for(k in 1:dim(driver.tab)[1]){
    ind<- which(grepl(driver.tab$POS[k], mut.vec, fixed = TRUE))
    mut.vec.driver[ind]<-1
  }
  
  #-- Spot-level driver status
  if(sum(mut.vec.driver)==1){
    m<-dat.mut[, which(mut.vec.driver==1)]
    df_loc$Driver<-(m==1)*1
  } else if (sum(mut.vec.driver)>1){
    m<-dat.mut[, which(mut.vec.driver==1)]
    df_loc$Driver<-(rowSums(m==1)*1)
  } else if(sum(mut.vec.driver)==0){
    df_loc$Driver=0
  }
  
  l<-dim(dat.mut)[2]
  df_loc[,6:(5+l)]<-dat.mut
  df<-rbind(df, df_loc)
  
}



#-----------------------------
#-- PLOTTING PREP
#-----------------------------

#-- Palette for Mutation calls
col_fun1 = colorRamp2(c(-1, 0, 0.5, 1), c("#EFEFEF", pal1[1], pal1[3], pal1[2]))

#-- Palette for Driver mutation count
col_fun2 = colorRamp2(c(0,max(df$Driver)), c("white", pal1[10] ))

#-- Palette for histology
col_fun3 = colorRamp2(c(0, 1, 2, 3), c(pal1[7], pal1[5], pal1[6], pal1[9]))

mat3<-rep(NA, dim(df)[1])
mat3[which(df$Histology== "Normal"  )]<-0
mat3[which(df$Histology== "Non-Atypical" | df$Histology==  "Atypical"   )]<-1
mat3[which(df$Histology== "DCIS"  )]<-2
mat3[which(df$Histology== "Invasive"  )]<-3

#-- Palette for tumor type
col_fun4 = colorRamp2(c(0, 1), c(pal1[3],  pal1[9]), transparency = .6)

mat4<-rep(NA, dim(df)[1])
mat4[which(df$Type=="Pure")]<-0
mat4[which(df$Type=="Mixed")]<-1

mat1<-as.matrix(df[,6:105])

ll<-which.max(colSums(is.na(mat1)))

mat1<-mat1[,1:(ll-1)]

mat1[is.na(mat1)]<--1

split<-as.vector(df$Tumor)
split<- sub(".*_", "", split)
split <- sub("\\.", "", split)


### Plot

circos.clear()


circlize_plot = function() {
  
  #-- Initialize
  circos.heatmap.initialize(mat1, split)
  
  #-- Mutations
  circos.heatmap(mat1, split = split, col = col_fun1, track.height=.4 )
  
  
  #-- Numbers
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + convert_y(2, "mm"), 
                CELL_META$sector.index,
                facing = "bending.inside", cex = 0.8,
                adj = c(0.5, 0), niceFacing = TRUE)
  }, bg.border = NA)
  
  
  #-- Histology 
  circos.heatmap(mat3, split=split, col=col_fun3, track.height=.1)
  
  #-- Drivers
  circos.heatmap(df$Driver, split = split, col = col_fun2, track.height=.025)
  
  #-- Type
  circos.heatmap(mat4, split=split, col=col_fun4, track.height=.2)
  
  circos.clear()
}


#-- Save Figure 2B as a PDF
pdf("Figure_2B.pdf")
circlize_plot()
dev.off()



