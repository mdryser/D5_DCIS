#-----------------------------------------------
#-----------------------------------------------
#-- Spot-specific MAF distributions
#-- Use for paper Figure 3A
#-- Updated: 09/13/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------


#-----------------------------------------------
#-- Load packages
#-----------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)

#-----------------------------------------------
#-- Initializing
#-----------------------------------------------

#-- Public mutation definition threshold
public.thr<-.9
#-- Mutation present threshold
post.mut.thr<-.95
#-- Fraction of spots covered to assess public mutation
public.spot<-.7

#-- The tumor to assay
tum="Tumor_173."

general_name_Mutations     <- "Data/Cleaned_Mutation_calls/"
general_name_MAF   <- "Data/Cleaned_VAF/"

#-- Import the files
Mut_files   <- list.files(path=general_name_Mutations, full.names = FALSE)
MAF_files  <- list.files(path=general_name_MAF, full.names = FALSE)

#-- The tumors
tumors<- sub(".*_", "", Mut_files)
tumors<- sub("csv.*", "", tumors)
tumors<-paste("Tumor_", tumors, sep="")

#-- Tumor-specific summaries
master.sum<-data.frame(tumor=tumors,
                       spots=rep(0,length(tumors)),
                       mut.total=rep(0,length(tumors)),
                       mut.public=rep(0,length(tumors)),
                       mut.private=rep(0,length(tumors)),
                       mut.unk=rep(0,length(tumors)),
                       mut.reads=rep(0,length(tumors))
)

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

#-- Identify clonal and subclonal mutations

read<-as.vector(colSums(!is.na(dat.mut))) # number of valid reads in that mutation

called<-as.vector(colSums(dat.mut>=post.mut.thr, na.rm = TRUE)) # number of spots the mutation is present

public.ind<-which(((called/read)>=public.thr) & (read>=public.spot*dim(dat.mut)[1]))
private.ind<-which(((called/read)<public.thr) & (read>=public.spot*dim(dat.mut)[1]))
NA.ind<-which(read<public.spot*dim(dat.mut)[1])

#----------------------------------------
#-- Select spots and plot the MAFs
#----------------------------------------

s1<-10
s2<-15
s3<-20
s4<-30


spot.1<-data.frame(spot=rep(rownames(MAF.redox)[s1] , dim(MAF.redox)[2]),
                   MAF=MAF.redox[s1,], mutation=colnames(MAF.redox))

spot.2<-data.frame(spot=rep(rownames(MAF.redox)[s2] , dim(MAF.redox)[2]),
                   MAF=MAF.redox[s2,], mutation=colnames(MAF.redox))

spot.3<-data.frame(spot=rep(rownames(MAF.redox)[s3] , dim(MAF.redox)[2]),
                   MAF=MAF.redox[s3,], mutation=colnames(MAF.redox))

spot.4<-data.frame(spot=rep(rownames(MAF.redox)[s4] , dim(MAF.redox)[2]),
                   MAF=MAF.redox[s4,], mutation=colnames(MAF.redox))

out<-rbind(spot.1, spot.2, spot.3, spot.4)

mut1<-which(out$mutation=="chr5.174948896.SFXN1")
mut2<-which(out$mutation=="chr4.164061514.NAF1")
mut3<-which(out$mutation!="chr4.164061514.NAF1" & out$mutation!="chr5.174948896.SFXN1")

out$mutation[mut1]<-1
out$mutation[mut2]<-2
out$mutation[mut3]<-3

cols <- c("1" = "#0073C2FF" , "2" = "#CD534CFF", "3" = "#868686FF")

out<- out %>% 
  filter(!is.na(MAF))

ggplot(out, aes(x = spot, y = MAF, colour=mutation)) + 
  geom_beeswarm(cex=3, size=4, show.legend = FALSE)+
  ylim(0,1)+
  xlab("Spot ID (DCIS-173)") +
  scale_colour_manual(values=cols) + 
  ylab("Variant allele frequency") +
  theme_bw()+
  theme(text = element_text(size = 20))  

#----------------------------------------
#-- Save Figure 3A to PDF
#----------------------------------------

ggsave(filename = "Figure_3A.pdf",
       width=20,
       height=20,
       units="cm")  

