#--------------------------------------------------------------------
#--------------------------------------------------------------------
#-- Spatial Spot Maps
#--------------------------------------------------------------------
#--------------------------------------------------------------------
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggbeeswarm)
library(ggsci)
#color palette
pal1 <- pal_jco("default")(10)


# The tumor to assay
tum="Tumor_53."


# Scale factors
f=3.5 # z-axis
f1=1.2 # y-axis
f2=1.2 # x-axis
f3=1.2 # slide label
f4=-.6 #spot label adjustment
bar.x=-1.1 # left of scale bar
bar.y=-.7

#--------------------------------------------------------------------------------------------------------------------------------------------
#---------------------------- Import the relevant files
#--------------------------------------------------------------------------------------------------------------------------------------------
#-- Load the paths
general_name_Mutations     <- "Data/Cleaned_Mutation_calls/"
general_name_MAF   <- "Data/Cleaned_VAF/"

Mut_files   <- list.files(path=general_name_Mutations, full.names = FALSE)


#-- Fetch the  data
ind<-grep(tum ,Mut_files, fixed=TRUE, value=TRUE)
dat.pos<-read.csv(file=paste(general_name_Mutations,ind,sep="/")) %>%
  filter(Spot!="Control") %>%
  remove_rownames() %>% 
  select(Spot, x.pos, y.pos, z.pos, Histology, Slide)



# The slides
u<-unique(dat.pos$Slide)

# Start-stop in x-axis direction
ind<-which(dat.pos$Slide==u[1])
max.1<-max((dat.pos$x.pos[ind]))
min.1<-min((dat.pos$x.pos[ind]))
# Z-axis position
z.1<-dat.pos$z.pos[ind[1]]*f


# Start-stop in x-axis direction
ind<-which(dat.pos$Slide==u[2])
max.2<-max((dat.pos$x.pos[ind]))
min.2<-min((dat.pos$x.pos[ind]))
# Z-axis position
z.2<-dat.pos$z.pos[ind[1]]*f

# Start-stop in x-axis direction
ind<-which(dat.pos$Slide==u[3])
max.3<-max((dat.pos$x.pos[ind]))
min.3<-min((dat.pos$x.pos[ind]))
# Z-axis position
z.3<-dat.pos$z.pos[ind[1]]*f


max.x<-max(abs(dat.pos$x.pos))*f2
max.y<-max(abs(dat.pos$y.pos))*f1


ggplot(dat.pos, aes(x=x.pos+f*z.pos, y=y.pos, label=Spot))+
  annotate("rect", xmin = z.1-max.x, xmax = z.1+max.x, ymin = -max.y, ymax = max.y,
           alpha = .3,fill = pal1[3]) +
  annotate("text", x = z.1, y = f3*max.y, label = u[1], size=3) +
  annotate("rect", xmin = z.2-max.x, xmax = z.2+max.x, ymin = -max.y, ymax = max.y,
           alpha = .3,fill =  pal1[3])+
  annotate("text", x = z.2, y = f3*max.y, label = u[2], size=3) +
  annotate("rect", xmin = z.3-max.x, xmax = z.3+max.x, ymin = -max.y, ymax = max.y,
           alpha = .3,fill =  pal1[3])+
  annotate("text", x = z.3, y = f3*max.y, label = u[3], size=3) +
  geom_point(aes(colour=factor(Histology)), size=2 )  +
  geom_text(hjust=f4, vjust=0, size=2, color="white") +
  scale_colour_manual(values=c("Normal"=pal1[7], "Benign"=pal1[5], "DCIS"=pal1[6], "Invasive"=pal1[9]))+
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
  geom_segment(aes(x=bar.x, y=bar.y, xend=bar.x+.5, yend=bar.y))+
  annotate(geom="text", x=bar.x+.2, y=bar.y+.2, label="5mm",
           color="black", size=2)



ggsave(paste(tum, "Spatial.Map.pdf", sep=""),
       path="Plots/",
       width=40,
       height=20,
       units="cm")
