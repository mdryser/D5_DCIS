#-----------------------------------------------
#-----------------------------------------------
#-- Mutation map for DCIS-168
#-- Use for paper Figure 4C
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
library(ggsci)


#-----------------------------------------------
#-- Initials
#-----------------------------------------------

#-- Color palette
pal1 <- pal_jco("default")(10)

#-- Mutation present threshold
post.mut.thr<-.95

#-- The tumor to assay
tum="Tumor_168."

general_name_Mutations     <- "Data/Cleaned_Mutation_calls/"
general_name_MAF   <- "Data/Cleaned_VAF/"

#-- Import relevant files
Mut_files   <- list.files(path=general_name_Mutations, full.names = FALSE)
MAF_files  <- list.files(path=general_name_MAF, full.names = FALSE)


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


ind<-grep(tum ,Mut_files, fixed=TRUE, value=TRUE)
dat.pos<-read.csv(file=paste(general_name_Mutations,ind,sep="/")) %>%
  filter(Spot!="Control") %>%
  filter(Histology=="DCIS") %>%
  remove_rownames() %>% 
  column_to_rownames(var = "Spot") 

#-----------------------------------------------
#-- Choices
#-----------------------------------------------

#-- Mutation choice
#-- For 24, use 17 for scattered, and 7 for contiguous 
#-- For 91, use 33 (chr2.238484162.RAB17) for scattered
#-- For 168, use 26 ("chr4...1836018.LETM1") for contiguous
mut.ind<-26


#-- Scaling of the z-axis
f=3.5

#-- Slide scalings
u<-unique(dat.pos$Slide)

ind<-which(dat.pos$Slide==u[1])
max.1<-max((dat.pos$x.pos[ind]))
min.1<-min((dat.pos$x.pos[ind]))
z.1<-dat.pos$z.pos[ind[1]]*f

ind<-which(dat.pos$Slide==u[2])
max.2<-max((dat.pos$x.pos[ind]))
min.2<-min((dat.pos$x.pos[ind]))
z.2<-dat.pos$z.pos[ind[2]]*f

ind<-which(dat.pos$Slide==u[3])
max.3<-max((dat.pos$x.pos[ind]))
min.3<-min((dat.pos$x.pos[ind]))
z.3<-dat.pos$z.pos[ind[3]]*f

maxx<-max(abs(dat.pos$x.pos))*1.2
max.y<-max(abs(dat.pos$y.pos))*2

#-- The mutation
ind.yes<-which(dat.mut[,mut.ind]>post.mut.thr)
ind.no<-which(dat.mut[,mut.ind]<(1-post.mut.thr))

dat.pos$mut<-"NA"
dat.pos$mut[ind.yes]<- "yes"
dat.pos$mut[ind.no]<- "no"


#-- The plot
ggplot(dat.pos, aes(x=x.pos+f*z.pos, y=y.pos))+
  annotate("rect", xmin = z.1-maxx, xmax = z.1+maxx, ymin = -max.y, ymax = max.y,
           alpha = .3,fill = pal1[3]) +
  annotate("text", x = z.1, y = 1.2*max.y, label = u[1]) +
  annotate("rect", xmin = z.2-maxx, xmax = z.2+maxx, ymin = -max.y, ymax = max.y,
           alpha = .3,fill =  pal1[3])+
  annotate("text", x = z.2, y = 1.2*max.y, label = u[2]) +
  annotate("rect", xmin = z.3-maxx, xmax = z.3+maxx, ymin = -max.y, ymax = max.y,
           alpha = .3,fill =  pal1[3])+
  annotate("text", x = z.3, y = 1.2*max.y, label = u[3]) +
  geom_point(aes(colour=factor(mut)), size=1.8 )  +
  scale_colour_manual(values=c("yes"=pal1[2], "no"=pal1[1], "NA"=pal1[3]))+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  coord_fixed(ratio=1)+  
  scale_x_continuous("Long axis (cm)", breaks=c( z.1,z.2,z.3), labels=c( c(z.1,z.2,z.3)/f))+
  ylab("y-axis (cm)")+
  theme(axis.line.y=element_line(color='white'), 
        axis.text.y = element_blank() ,
        axis.ticks.y = element_blank(),
        legend.position = "none")+
  ylab("") + 
  geom_segment(aes(x=0, y=-1.2, xend=.5, yend=-1.2))+
  annotate(geom="text", x=.25, y=-1.05, label="5mm",
           color="black", size=2)


#-- Save Figure 4C to PDF
ggsave("Figure_4C.pdf",
       width=20,
       height=10,
       units="cm")
