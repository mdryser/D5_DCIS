#-----------------------------------------------
#-----------------------------------------------
#-- Mutation Venn Diagrams
#-- Use for paper Figure 6B, 6E
#-- Updated: 09/13/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------

#-----------------------------------------------
#-- Load packages
#-----------------------------------------------
library(readr)
library(ggplot2)
library(dplyr)
library(Rtsne)
library(ggsci)
library(VennDiagram)
library(ggVennDiagram)
library(Rcpp)


#-----------------------------------------------
#-- Preparations
#-----------------------------------------------

#-- Posterior probability cutoff for calling a mutation present
prob.mut.thr<-.95

#-- Color palette
pal1 <- pal_jco("default")(10)

#-- Not including #222 because there are only 2 DCIS spots
Tumornames <- c("03","17","24","53","65","66","88","91","168","173","191","211","213","221","222","225","227","286","315")


#-----------------------------------------------
#-- Tumor Processing
#------------------------------------------------
f.proc<-function(DCIS_name, prob.mut.thr){
  
  #- Minor pre-processing
  DCIS_path <- paste("Data/Cleaned_Mutation_calls/CleanedBayesian_Tumor_", DCIS_name, ".csv", sep = "")
  Tumor <- read.csv(DCIS_path)
  Tumor$X <- NULL
  colnames(Tumor) <-  gsub(" ", "", colnames(Tumor), fixed = TRUE)
  
  #-- Adding Histology_3cat: a pre-DCIS histology category
  Tumor <- Tumor %>% mutate(Histology_3cat= 
                              ifelse(Histology %in% c("Atypical", "Non-Atypical","Normal"), 
                                     "preDCIS", Histology))
  
  #-- Adding Histology_4cat: normal-benign-dcis-invasiv
  Tumor <- Tumor %>% mutate(Histology_4cat= 
                              ifelse(Histology %in% c("Atypical", "Non-Atypical"), 
                                     "Benign", Histology))
  
  #-- Only select the columns with mutations
  nam<-grep('chr', colnames(Tumor), value=TRUE) 
  
  #--  Identify mutations where we know for sure they are present
  Tumor[, nam] <- ifelse(Tumor[, nam]>prob.mut.thr, 1, 0)
  
  #- Return the tumor
  return(list(Tumor, nam))
  
}



#-----------------------------------------------
#--- DCIS-03: Normal-DCIS
#------------------------------------------------
DCIS_name<-"03"
out<-f.proc(DCIS_name, prob.mut.thr)
Tumor<-out[[1]]
nam<-out[[2]]

Mut.sum.2cat <- list(preDCIS=0,
                     DCIS=0)

#-- Find mutations in each category
ind<-which(Tumor$Histology_3cat=="preDCIS")
Mut.sum.2cat$preDCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

ind<-which(Tumor$Histology_3cat=="DCIS")
Mut.sum.2cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

#-- The Venn diagram
ggVennDiagram(Mut.sum.2cat,
              edge_lty = 0,
              set_size = 8,
              label="count",
              label_size=15,
              category.names = c("Normal","DCIS"))+
  scale_fill_gradient(low="white",high = pal1[4])+
  theme(legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=16))+ #change legend text font size
  theme(text = element_text(size = 16))

ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
       width=12,
       height=12,
       units="cm")

a <- lengths(Mut.sum.2cat)
print(a)


#-----------------------------------------------
#--- DCIS-17: Normal-DCIS
#------------------------------------------------
DCIS_name<-"17"
out<-f.proc(DCIS_name, prob.mut.thr)
Tumor<-out[[1]]
nam<-out[[2]]

Mut.sum.2cat <- list(preDCIS=0,
                     DCIS=0)

#-- Find mutations in each category
ind<-which(Tumor$Histology_3cat=="preDCIS")
Mut.sum.2cat$preDCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

ind<-which(Tumor$Histology_3cat=="DCIS")
Mut.sum.2cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

#-- The Venn diagram
ggVennDiagram(Mut.sum.2cat,
              edge_lty = 0,
              set_size = 8,
              label="count",
              label_size=15,
              category.names = c("Normal","DCIS"))+
  scale_fill_gradient(low="white",high = pal1[4])+
  theme(legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=16))+ #change legend text font size
  theme(text = element_text(size = 16))

ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
       width=12,
       height=12,
       units="cm")

a <- lengths(Mut.sum.2cat)
print(a)


#-----------------------------------------------
#--- DCIS-24: Benign-DCIS
#------------------------------------------------
DCIS_name<-"24"
out<-f.proc(DCIS_name, prob.mut.thr)
Tumor<-out[[1]]
nam<-out[[2]]

Mut.sum.2cat <- list(preDCIS=0,
                     DCIS=0)

#-- Find mutations in each category
ind<-which(Tumor$Histology_3cat=="preDCIS")
Mut.sum.2cat$preDCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

ind<-which(Tumor$Histology_3cat=="DCIS")
Mut.sum.2cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

#-- The Venn diagram
ggVennDiagram(Mut.sum.2cat,
              edge_lty = 0,
              set_size = 8,
              label="count",
              label_size=15,
              category.names = c("Benign", "DCIS"))+
  scale_fill_gradient(low="white",high = pal1[4])+
  theme(legend.key.size = unit(2, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=16))+ #change legend text font size
  theme(text = element_text(size = 16))

ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
       width=12,
       height=12,
       units="cm")

a <- lengths(Mut.sum.2cat)
print(a)

#-----------------------------------------------
#--- DCIS-213: Normal-Benign-DCIS
#-----------------------------------------------
DCIS_name<-"213"
out<-f.proc(DCIS_name, prob.mut.thr)
Tumor<-out[[1]]
nam<-out[[2]]

Tumor <- Tumor %>% filter(Spot != 6)

Mut.sum.3cat <- list(Benign=0,
                     Normal=0,
                     DCIS=0)

#-- Find mutations in each category
ind<-which(Tumor$Histology_4cat=="Normal")
Mut.sum.3cat$Normal<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

ind<-which(Tumor$Histology_4cat=="Benign")
Mut.sum.3cat$Benign<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

ind<-which(Tumor$Histology_4cat=="DCIS")
Mut.sum.3cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]


#-- The Venn diagram
ggVennDiagram(Mut.sum.3cat,
              edge_lty = 0,
              set_size = 8,
              label="count",
              label_size=15,
              category.names = c( "Benign (atypical)","Normal", "DCIS")) + 
  scale_fill_gradient(low="white",high = pal1[4])+
  theme(legend.key.size = unit(4, 'cm'), #change legend key size
        legend.key.height = unit(2, 'cm'), #change legend key height
        legend.key.width = unit(2, 'cm'), #change legend key width
        legend.title = element_text(size=20), #change legend title font size
        legend.text = element_text(size=30))+ #change legend text font size
  theme(text = element_text(size = 30))



ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
       width=20,
       height=20,
       units="cm")

a <- lengths(Mut.sum.3cat)
print(a)




#---------------------------------------------
#--- DCIS 222: Benign(-)-Benign(+)-DCIS
#---------------------------------------------
DCIS_name<-"222"
out<-f.proc(DCIS_name, prob.mut.thr)
Tumor<-out[[1]]
nam<-out[[2]]

  Mut.sum.3cat <- list(Atypia=0,
                       NonAtypia=0,
                       DCIS=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology=="Non-Atypical")
  Mut.sum.3cat$NonAtypia<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="Atypical")
  Mut.sum.3cat$Atypia<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="DCIS")
  Mut.sum.3cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.3cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=15,
                category.names = c( "Benign (w atypia)","Benign", "DCIS")) + 
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(4, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(2, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=30))+ #change legend text font size
    theme(text = element_text(size = 30))
  
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=20,
         height=20,
         units="cm")
  
  a <- lengths(Mut.sum.3cat)
  print(a)
  


  
  #-----------------------------------------------
  #--- DCIS-315: Benign-DCIS
  #------------------------------------------------
  DCIS_name<-"315"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  
  Mut.sum.2cat <- list(preDCIS=0,
                       DCIS=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology_3cat=="preDCIS")
  Mut.sum.2cat$preDCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology_3cat=="DCIS")
  Mut.sum.2cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.2cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=15,
                category.names = c("Benign","DCIS"))+
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=16), #change legend title font size
          legend.text = element_text(size=16))+ #change legend text font size
    theme(text = element_text(size = 16))
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=12,
         height=12,
         units="cm")
  
  a <- lengths(Mut.sum.2cat)
  print(a)
  
  
  #-----------------------------------------------
  #--- DCIS-53: Benign (w and without -DCIS-Invasive
  #------------------------------------------------
  DCIS_name<-"53"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  Mut.sum.4cat <-list(Normal=0,
                      Benign=0,
                      DCIS=0,
                      Invasive=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology=="Non-Atypical")
  Mut.sum.4cat$Normal<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="Atypical")
  Mut.sum.4cat$Benign<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="DCIS")
  Mut.sum.4cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="Invasive")
  Mut.sum.4cat$Invasive<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.4cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=8,
                category.names = c("Benign (-)","Benign (+)","DCIS", "Invasive" )) + 
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(4, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(2, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=30)) #change legend text font size
  
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=35,
         height=25,
         units="cm")
  
  a <- lengths(Mut.sum.4cat)
  print(a)
  
  
  
  #-----------------------------------------------
  #--- DCIS-65: Normal-DCIS-Invasive
  #------------------------------------------------
  DCIS_name<-"65"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  Mut.sum.3cat <- list(DCIS=0,
                       Normal=0,
                       Invasive=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology=="Normal")
  Mut.sum.3cat$Normal<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="DCIS")
  Mut.sum.3cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="Invasive")
  Mut.sum.3cat$Invasive<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.3cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=15,
                category.names = c( "DCIS","Normal", "Invasive")) + 
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(4, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(2, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=30))+ #change legend text font size
    theme(text = element_text(size = 30))
  
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=20,
         height=20,
         units="cm")
  
  a <- lengths(Mut.sum.3cat)
  print(a)
  
  
  
  #-----------------------------------------------
  #--- DCIS-66: Normal-Benign-DCIS-Invasive
  #------------------------------------------------
  DCIS_name<-"66"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  Mut.sum.4cat <-list(Normal=0,
                      Benign=0,
                      DCIS=0,
                      Invasive=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology=="Normal")
  Mut.sum.4cat$Normal<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology_4cat=="Benign")
  Mut.sum.4cat$Benign<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="DCIS")
  Mut.sum.4cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="Invasive")
  Mut.sum.4cat$Invasive<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.4cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=8,
                category.names = c("Normal","Benign","DCIS", "Invasive" )) + 
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(4, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(2, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=30)) #change legend text font size
  
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=35,
         height=25,
         units="cm")
  
  
  a <- lengths(Mut.sum.4cat)
  print(a)
  
  
  
  #-----------------------------------------------
  #--- DCIS-88: DCIS-Invasive
  #------------------------------------------------
  DCIS_name<-"88"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  
  Mut.sum.2cat <- list(DCIS=0,
                       Invasive=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology_3cat=="DCIS")
  Mut.sum.2cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology_3cat=="Invasive")
  Mut.sum.2cat$Invasive<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.2cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=15,
                category.names = c("DCIS", "Invasive"))+
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=16), #change legend title font size
          legend.text = element_text(size=16))+ #change legend text font size
    theme(text = element_text(size = 16))
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=12,
         height=12,
         units="cm")
  
  a <- lengths(Mut.sum.2cat)
  print(a)
  
  
  
  #-----------------------------------------------
  #--- DCIS-91: DCIS-Invasive
  #------------------------------------------------
  DCIS_name<-"91"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  Mut.sum.2cat <- list(DCIS=0,
                       Invasive=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology_3cat=="DCIS")
  Mut.sum.2cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology_3cat=="Invasive")
  Mut.sum.2cat$Invasive<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.2cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=15,
                category.names = c("DCIS", "Invasive"))+
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=16), #change legend title font size
          legend.text = element_text(size=16))+ #change legend text font size
    theme(text = element_text(size = 16))
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=12,
         height=12,
         units="cm")
  
  a <- lengths(Mut.sum.2cat)
  print(a)
  
  
  
  #-----------------------------------------------
  #--- DCIS-168: DCIS-Invasive
  #------------------------------------------------
  DCIS_name<-"168"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  Mut.sum.2cat <- list(DCIS=0,
                       Invasive=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology_3cat=="DCIS")
  Mut.sum.2cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology_3cat=="Invasive")
  Mut.sum.2cat$Invasive<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.2cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=15,
                category.names = c("DCIS", "Invasive"))+
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=16), #change legend title font size
          legend.text = element_text(size=16))+ #change legend text font size
    theme(text = element_text(size = 16))
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=12,
         height=12,
         units="cm")
  
  a <- lengths(Mut.sum.2cat)
  print(a)
  
  

  #-----------------------------------------------
  #--- DCIS-173: Normal-DCIS-Invasive
  #------------------------------------------------
  DCIS_name<-"173"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  Mut.sum.3cat <- list(DCIS=0,
                       Normal=0,
                       Invasive=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology=="Normal")
  Mut.sum.3cat$Normal<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="DCIS")
  Mut.sum.3cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology=="Invasive")
  Mut.sum.3cat$Invasive<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]

  #-- The Venn diagram
  ggVennDiagram(Mut.sum.3cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=15,
                category.names = c( "DCIS","Normal", "Invasive")) + 
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(4, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(2, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=30))+ #change legend text font size
    theme(text = element_text(size = 30))
  
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=20,
         height=20,
         units="cm")
  
  a <- lengths(Mut.sum.3cat)
  print(a)
  
  
  
  
  #-----------------------------------------------
  #--- DCIS-191: Benign-DCIS
  #------------------------------------------------
  DCIS_name<-"191"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  
  Mut.sum.2cat <- list(Benign=0,
                       DCIS=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology=="Non-Atypical")
  Mut.sum.2cat$Benign<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology_3cat=="DCIS")
  Mut.sum.2cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.2cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=15,
                category.names = c("Benign", "DCIS"))+
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(2, 'cm'), #change legend key size
          legend.key.height = unit(1, 'cm'), #change legend key height
          legend.key.width = unit(1, 'cm'), #change legend key width
          legend.title = element_text(size=16), #change legend title font size
          legend.text = element_text(size=16))+ #change legend text font size
    theme(text = element_text(size = 16))
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=12,
         height=12,
         units="cm")
  
  a <- lengths(Mut.sum.2cat)
  print(a)
  
  
  
  
  #-----------------------------------------------
  #--- DCIS-286: Normal-Benign-DCIS-Invasive
  #------------------------------------------------
  DCIS_name<-"286"
  out<-f.proc(DCIS_name, prob.mut.thr)
  Tumor<-out[[1]]
  nam<-out[[2]]
  
  # remove the normal problem spot
  
  Tumor <- Tumor %>% filter(Histology!="Normal")
  
  Mut.sum.3cat <- list(DCIS=0,
                       Benign=0,
                       Invasive=0)
  
  #-- Find mutations in each category
  ind<-which(Tumor$Histology_4cat=="Benign")
  Mut.sum.3cat$Benign<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology_4cat=="DCIS")
  Mut.sum.3cat$DCIS<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  ind<-which(Tumor$Histology_4cat=="Invasive")
  Mut.sum.3cat$Invasive<- nam[colSums(Tumor[ind,nam], na.rm = TRUE)>0]
  
  #-- The Venn diagram
  ggVennDiagram(Mut.sum.3cat,
                edge_lty = 0,
                set_size = 8,
                label="count",
                label_size=8,
                category.names = c("DCIS","Benign", "Invasive" )) + 
    scale_fill_gradient(low="white",high = pal1[4])+
    theme(legend.key.size = unit(4, 'cm'), #change legend key size
          legend.key.height = unit(2, 'cm'), #change legend key height
          legend.key.width = unit(2, 'cm'), #change legend key width
          legend.title = element_text(size=20), #change legend title font size
          legend.text = element_text(size=30)) #change legend text font size
  
  
  ggsave(file=paste("Venn_diagrams/DCIS_", DCIS_name, ".pdf", sep=""),
         width=35,
         height=25,
         units="cm")
  
  
  a <- lengths(Mut.sum.3cat)
  print(a)
  
  
