#-----------------------------------------------
#-----------------------------------------------
#-- Mutation diameters in DCIS portions
#-- Use for paper Figure 4A
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
library(stringr)
library(ggsci)

#-----------------------------------------------
#-- Initials
#-----------------------------------------------

#-- Color palette
pal1 <- pal_jco("default")(10)

#-- Public mutation definition threshold
public.thr<-.9
#-- Mutation present threshold
post.mut.thr<-.95
#-- Fraction of spots covered to assess public mutation
public.spot<-.75


#-- Load the paths
general_name_Mutations     <- "Data/Cleaned_Mutation_calls/"
general_name_MAF   <- "Data/Cleaned_VAF/"


#-----------------------------------------------
#-- Auxiliary diameter function
#-----------------------------------------------
# INPUT: mut.vec is the posterior probability of a mutation being present in each spot
#        pos.vec: (x,y,z)-coordinates of the spots
#        post.mut.thr: used to call mutations
# OUTPUT: diam=diameter of the mutation
#         diam.max=maximum diameter for this mutation

f.diam <- function(mut.vec,pos.vec, post.mut.thr) {
  
  #-- Mutation calls
  out<-mut.vec*NA
  out[mut.vec>post.mut.thr]<-1
  out[mut.vec<(1-post.mut.thr)]<-0
  
  #-- Number of spots that contain the mutation
  ind<-which(out==1)
  
  #-- Record mutation diameter if >0 spots contain the mutation
  if(length(ind)>1){
    diam= max(dist(pos.vec[ind,]))
  } else if(length(ind)==1){
    diam=0
  
  #-- Record NA if no spot contains the mutation
  } else {
    diam=NA
  }
  
  #-- Diameter of spots that have a valid mutation read
  ind<-which(!is.na(out))
  diam.max= max(dist(pos.vec[ind,]))
  
  return(c(diam, diam.max))
}


#-----------------------------------------------
#-- Auxiliary function to determine public vs private (aka restricted)
#-----------------------------------------------
# INPUT: mut.vec is the posterior probability of a mutation being present in each spot
#        pos.vec: (x,y,z)-coordinates of the spots
#        post.mut.thr: used to call mutations
#        public.thr: minimum fractions covered to be classified as public
# OUTPUT: label: public vs  private

f.type<- function(mut.vec, post.mut.thr, public.thr, public.spot) {
  
  #-- Call the mutation
  out<-mut.vec*NA
  out[mut.vec>post.mut.thr]<-1
  out[mut.vec<(1-post.mut.thr)]<-0
  
  read<-sum(!is.na(out))
  called<-sum(out==1, na.rm = TRUE)
  
  if(((read/length(out))>=public.spot)& (called>0)){
    public_private=ifelse(called/read>public.thr, "public", "private")
  } else {
    public_private=NA
  }
  return(public_private)
}

#-----------------------------------------------
#-- Auxiliary Spot function 
#-----------------------------------------------
# INPUT: mut.vec is the posterior probability of a mutation being present in each spot
#        post.mut.thr: used to call mutations
# OUTPUT: spots: number of spots occupied by the mutation
#         spots.max: max number of spots that could be occupied (only valid spots for this mutation)

f.spots <- function(mut.vec, post.mut.thr) {
  
  out<-mut.vec*NA
  out[mut.vec>post.mut.thr]<-1
  out[mut.vec<(1-post.mut.thr)]<-0
  
  #-- Number of spots with mutation
  spots<-sum(out==1, na.rm = TRUE)
  
  #-- Max possible number of spots with mutation
  spots.max<-sum(!is.na(out))
  
  return(c(spots, spots.max))
}

#-----------------------------------------------
#-- Import the relevant files
#-----------------------------------------------

Mut_files   <- list.files(path=general_name_Mutations, full.names = FALSE)
MAF_files  <- list.files(path=general_name_MAF, full.names = FALSE)


tumors<- sub(".*_", "", Mut_files)
tumors<- sub("csv.*", "", tumors)
tumors<-paste("Tumor_", tumors, sep="")

ind<-which(tumors=="Tumor_221.")
tumors<-tumors[-ind]
# Compiling all the MAFs
master.diam<-data.frame(tumor="a",
                        name="b",
                        diam=0,
                        diam.max=0,
                        diam.global=0,
                        spots=0,
                        spots.max=0,
                        type="p")


#-----------------------------------------------
#-- Loop over all tumors
#-----------------------------------------------

for(tum in tumors ){
  
  #-- Fetch and clean the  mutation data
  ind<-grep(tum ,Mut_files, fixed=TRUE, value=TRUE)
  dat.mut<-read.csv(file=paste(general_name_Mutations,ind,sep="/"),check.names=FALSE)
  
  colnames(dat.mut)<-str_replace_all(colnames(dat.mut), fixed("-"), ".")
  colnames(dat.mut)<-str_replace_all(colnames(dat.mut), fixed(" "), "")
  
 dat.mut <- dat.mut %>%  select(-1) %>%
    filter(Spot!="Control") %>%
    filter(Histology=="DCIS") 
  
  #-- Create a streamlined mutation only data frame
  dat.mut.only<-dat.mut %>% remove_rownames() %>% 
    column_to_rownames(var = "Spot") %>%
    select(contains("chr"))
  
  #-- Calculate the maximum diameter in the tumor
  diam.global<- max(dist(dat.mut[,c("x.pos","y.pos","z.pos")]))
  
  #-- Calculate the diameter for each mutation
  diam.sum<-apply(dat.mut.only, 2, function(x) f.diam(x,dat.mut[,c("x.pos","y.pos","z.pos")], post.mut.thr))
  
  #-- Determine whether a mutation is public or restricted
  type<-apply(dat.mut.only, 2, function(x) f.type(x, post.mut.thr, public.thr, public.spot))
  
  # Determine the number of spots visited by the mutation
  spots<-apply(dat.mut.only, 2, function(x) f.spots(x, post.mut.thr))
  
  local.diam<-data.frame(tumor=rep(tum, dim(diam.sum)[2]),
                         name=names(diam.sum[1,]),
                         diam=diam.sum[1,],
                         diam.max=diam.sum[2,],
                         diam.global=rep(diam.global,dim(diam.sum)[2]),
                         spots=spots[1,],
                         spots.max=spots[2,],
                         type=type)
  
  master.diam<-rbind(master.diam, local.diam)
  
}

master.diam <- master.diam[-1,]


#-----------------------------------------------
#-- Driver mutation annotation
#-----------------------------------------------

#-- Convert master.maf to something more informative
master.diam$CHROM<-sub(".*chr", "", sub(":.*", "", master.diam$name))
master.diam$POS<-sub("(.*):.*", "\\1", master.diam$name)
master.diam$POS<-sub(".*:(.*)", "\\1", master.diam$POS)


MUT.annotation <- read.csv(file="Data/Annotations_SNV/SNV_annotations_complete.csv") %>%
  mutate(Driver=replace(Driver, Driver=="", "intergenic")) %>%
  rename(tumor=DCISsample)

tumor_rename<- function(old) {
  new<-paste("Tumor_", sub(".*d", "", old), ".", sep="") 
  return(new)
}

MUT.annotation$tumor<-sapply(MUT.annotation$tumor, tumor_rename)

table(MUT.annotation$Driver)

master.diam$Driver<-NA

for(k in 1: dim(master.diam)[1]){
  ind<-which(MUT.annotation$CHROM==master.diam$CHROM[k] &
               MUT.annotation$POS==master.diam$POS[k] &
               MUT.annotation$tumor==master.diam$tumor[k])
  
  if(length(ind)==0){print(k)}
  
  if(length(ind)>0){
    master.diam$Driver[k]<-MUT.annotation$Driver[ind]
  }
}

#-- Filter mutations that don't appear in any spot of the tumor
master.diam <- master.diam %>% 
  filter(!is.na(master.diam$diam)) %>% 
  mutate(diam.norm=diam/diam.max) %>%
  mutate(spots.norm=spots/spots.max)


#-----------------------------------------------
#-- Private/public and Driver/Passenger Statistics
#-----------------------------------------------
#-- Summaries
table(master.diam$Driver)
table(master.diam$type)

master.diam$Driver.bin<-ifelse(master.diam$Driver=="likely", "Driver", "Passenger")
table(master.diam$Driver.bin)

#-- Crass- tabulation
table(master.diam$Driver.bin, master.diam$type)

#-- Fisher exact test for difference
fisher.test(table(master.diam$Driver.bin, master.diam$type))


#-----------------------------------------------
#-- Diameter comparison: driver vs passenger
#-----------------------------------------------
ggplot(master.diam, aes(x = Driver.bin, y = diam.norm, fill = Driver.bin)) +
  geom_violin(alpha=1, scale="width") +
  scale_fill_manual(values = c("Driver" = pal1[4], "Passenger" = pal1[3])) +
  labs(y = "Normalized diameter (%)", x = "")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 15))+
  scale_x_discrete(labels= c(paste("Drivers (n=", sum(master.diam$Driver.bin=="Driver"), ")", sep=""),
                             c(paste("Passengers  (n=", sum(master.diam$Driver.bin=="Passenger"), ")", sep=""))))

ggsave(file="Figure_S6.pdf",
       width=10,
       height=10,
       units="cm")


driver<-master.diam$diam.norm[master.diam$Driver.bin=="Driver"]
passenger<-master.diam$diam.norm[master.diam$Driver.bin=="Passenger"]
t.test(driver, passenger)
wilcox.test(driver,passenger)


#-----------------------------------------------
## Fraction of spots covered: driver vs passenger
#-----------------------------------------------

ggplot(master.diam, aes(x = Driver.bin, y = spots.norm, fill = Driver.bin)) +
  geom_violin(alpha=.6, scale="width" ) +
  scale_fill_manual(values = c("Driver" = pal1[4], "Passenger" = pal1[3])) +
  labs(y = "Fraction of spots covered", x = "")+
  theme_bw()+
  theme(legend.position = "none")+
  theme(text = element_text(size = 15))+
  scale_x_discrete(labels= c(paste("Drivers (n=", sum(master.diam$Driver.bin=="Driver"), ")", sep=""),
                             c(paste("Passengers  (n=", sum(master.diam$Driver.bin=="Passenger"), ")", sep=""))))

ggsave(file="Figure_3C.pdf",
       width=10,
       height=10,
       units="cm")

driver<-master.diam$spots.norm[master.diam$Driver.bin=="Driver"]
passenger<-master.diam$spots.norm[master.diam$Driver.bin=="Passenger"]
median(driver)
median(passenger)
t.test(driver, passenger)
wilcox.test(driver,passenger)



#-----------------------------------------------
#-- Private mutation diameter figure (Fig 4A)
#-----------------------------------------------

master.diam.private <- master.diam  %>%
  filter(type=="private")

master.diam.private$tumor<-gsub(".*_","",master.diam.private$tumor)
master.diam.private$tumor<-str_replace(master.diam.private$tumor, "[.]", "")

a<- unique(master.diam.private[, c("tumor", "diam.global")])

#-- Order them 
ord<-c("03", "17", "24", "211", "213",  "225", "227", "315", "53", "65", "66", "88", "91", "168", "173", "191", "286")
a<-a[match(ord, a$tumor),]
a$tumor<-as.factor(a$tumor)
a$tumor <- factor(a$tumor, levels=ord)
a$type<-c(rep("Pure", 8), rep("Mixed",9))

#-- The Plot
ggplot() +
  geom_col(data=a,aes(tumor, diam.global, fill=type, alpha=.6))+
  geom_quasirandom(data=master.diam.private, aes(tumor, diam),varwidth = FALSE, size=1, groupOnX=TRUE, colour="black") +
  xlab("Tumor ID") + ylab("Mutation diameter (cm)") +
  theme_bw()+theme(text = element_text(size = 12))+
  scale_fill_manual(values=c(pal1[4], pal1[3]))

#-- Save Figure 4A to PDF
ggsave(file="Figure_4A.pdf",
       width=20,
       height=8,
       units="cm")





