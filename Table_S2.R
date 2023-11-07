#-----------------------------------------------
#-----------------------------------------------
#-- Mutational analyses across all tumors
#-- Use for Supplementary Table S2
#-- Updated: 09/13/2022 by mdr30
#-----------------------------------------------
#-----------------------------------------------

#----------------------------------------
#-- CHOICE: "All" SPOTS OR "DCIS" SPOTS?
#----------------------------------------
type<-"DCIS"


#----------------------------------
#-- Preparations
#----------------------------------

#-- Libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggsci)

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

#----------------------------------
#-- Import the relevant files
#----------------------------------

Mut_files   <- list.files(path=general_name_Mutations, full.names = FALSE)
MAF_files  <- list.files(path=general_name_MAF, full.names = FALSE)


tumors<- sub(".*_", "", Mut_files)
tumors<- sub("csv.*", "", tumors)
tumors<-paste("Tumor_", tumors, sep="")

# Compiling all the MAFs
master.maf<-data.frame(tumor="a",
                       Mutation="b",
                       MAF=0,
                       label="y")

# Tumor-specific summaries
master.sum<-data.frame(tumor=tumors,
                       spots=rep(0,length(tumors)),
                       mut.total=rep(0,length(tumors)),
                       mut.public=rep(0,length(tumors)),
                       mut.private=rep(0,length(tumors)),
                       mut.unk=rep(0,length(tumors)),
                       mut.reads=rep(0,length(tumors))
)


#----------------------------------
#-- MAIN LOOP OVER ALL TUMORS
#----------------------------------

for(tum in tumors ){
  
  #-- Fetch the mutation data
  ind<-grep(tum ,Mut_files, fixed=TRUE, value=TRUE)
  dat.mut<-read.csv(file=paste(general_name_Mutations,ind,sep="/"),check.names=FALSE)  %>%
    select(-1)
  
  #-- Clean the column names
  colnames(dat.mut)<-str_replace_all(colnames(dat.mut), fixed("-"), ".")
  colnames(dat.mut)<-str_replace_all(colnames(dat.mut), fixed(" "), "")
  
  
  #-- Filter out the DCIS spots if needed
  if(type=="All"){
    dat.mut = dat.mut %>%
      filter(Spot!="Control") %>%
      remove_rownames() %>% 
      column_to_rownames(var = "Spot") %>%
      select(contains("chr"))
  } else if(type=="DCIS"){
    dat.mut = dat.mut %>%
      filter(Spot!="Control") %>%
      filter(Histology=="DCIS") %>%
      remove_rownames() %>% 
      column_to_rownames(var = "Spot") %>%
      select(contains("chr"))
  } else {
    
    print('warning')
  }

  
  #-- Fetch MAF files
  ind<-grep(tum ,MAF_files, fixed=TRUE, value=TRUE)
  dat.MAF<-read.csv(file=paste(general_name_MAF, ind,sep="/"),check.names=FALSE) %>%
    select(-1) 
  
  #-- Clean the column names
  colnames(dat.MAF)<-str_replace_all(colnames(dat.MAF), fixed("-"), ".")
  colnames(dat.MAF)<-str_replace_all(colnames(dat.MAF), fixed(" "), "")
  
  
  #-- Filter out the DCIS spots if needed
  if(type=="All"){
  dat.MAF = dat.MAF %>%
    filter(Spot!="Control") %>%
    remove_rownames() %>% 
    column_to_rownames(var = "Spot") %>%
    select(contains("chr"))
  } else if(type=="DCIS"){
    dat.MAF = dat.MAF %>%
      filter(Spot!="Control") %>%
      filter(Histology=="DCIS") %>%
      remove_rownames() %>% 
      column_to_rownames(var = "Spot") %>%
      select(contains("chr"))
  } else {
    
    print('warning')
  }
    
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
  
  #-- Identify Public and Private mutations
  #-- Number of mutations that are found in each spot
  read<-as.vector(colSums(!is.na(dat.mut))) # number of valid reads in that mutation
  #-- Number of mutations present
  called<-as.vector(colSums(dat.mut>=post.mut.thr, na.rm = TRUE)) # number of spots the mutation is present
  
  #-- Record public mutations
  public.ind<-which(((called/read)>=public.thr) & (read>=public.spot*dim(dat.mut)[1]))
  #-- Record private mutations
  private.ind<-which(((called/read)<public.thr) & (read>=public.spot*dim(dat.mut)[1]))
  #-- Record NA mutations
  NA.ind<-which(read<public.spot*dim(dat.mut)[1])
  
  #-- Assemble the public mutations
  if(length(public.ind)>1){
    MAF.public <- as.data.frame(MAF.redox[,public.ind]) %>%
      pivot_longer(cols=everything(),names_to = "Mutation", values_to="MAF")%>%
      mutate(label="public") %>% 
      filter(!is.na(MAF))
  }
  
  #-- Special case of just one public mutation
  if(length(public.ind)==1){
    MAF.public <- data.frame(Mutation=rep(colnames(MAF.redox)[public.ind], dim(MAF.redox)[1] ), MAF=MAF.redox[,public.ind]) %>%
      mutate(label="public") %>% 
      filter(!is.na(MAF))
  }
  
  
  #-- Assemble the private mutations
  if(length(private.ind)>0){
    MAF.private <- as.data.frame(MAF.redox[,private.ind]) %>%
      pivot_longer(cols=everything(),names_to = "Mutation", values_to="MAF") %>%
      mutate(label="private") %>% 
      filter(!is.na(MAF))
  }
  
  
  #-- Assemble the NA mutations
  if(length(NA.ind)>0){
    MAF.NA <- as.data.frame(MAF.redox[,NA.ind]) %>%
      pivot_longer(cols=everything(),names_to = "Mutation", values_to="MAF") %>%
      mutate(label=NA) %>% 
      filter(!is.na(MAF))
  }
  
  if(length(NA.ind)==1){
    MAF.NA <- data.frame(Mutation=rep(colnames(MAF.redox)[NA.ind], dim(MAF.redox)[1] ), MAF=MAF.redox[,NA.ind]) %>%
      mutate(label=NA) %>% 
      filter(!is.na(MAF))
  }
  
  
  #-- Concatenate public and private mutations
  lpu=length(public.ind)
  lpr=length(private.ind)
  
  
  if(lpu==0 & lpr>0){
    MAF.out<-MAF.private
  } else if (lpu>0 & lpr==0){
    MAF.out<-MAF.public
  } else if (lpu>0 & lpr>0){
    MAF.out<-rbind(MAF.public, MAF.private)  
  }
  
  #-- Append the NAs
  if(length(NA.ind)>0){
    
    MAF.out<-rbind(MAF.out, MAF.NA)
  }
  
  MAF.out<- MAF.out %>%
    mutate(tumor=tum)
  
  master.maf<-rbind(master.maf,MAF.out)
  
  
  #-- Compile tumor summary sheet
  ind<-which(master.sum$tumor==tum)
  master.sum$mut.public[ind]<-lpu
  master.sum$mut.private[ind]<-lpr
  master.sum$mut.unk[ind]<-length(NA.ind)
  master.sum$mut.total[ind]<-length(NA.ind)+lpu+lpr
  master.sum$spots[ind]<-dim(dat.mut)[1]
  master.sum$mut.reads[ind]<-dim(MAF.out)[1]
  
  
}

#-- Remove the initializing row
master.maf<-master.maf[-1,]


#-------------------------------------
#-- Plot histogram private vs public (Figure S7A)
#-------------------------------------

ind<-which(is.na(master.maf$label))

ggplot(master.maf[-ind,], aes(x=MAF,  fill=label)) +
  scale_fill_manual(values=c("public"=pal1[3],"private"=pal1[1] ))+
  geom_histogram(alpha=0.6, position="identity", binwidth = .03)+
  coord_cartesian(xlim=c(0,1)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size = 18)) +
  ylab("Count")+
  xlab("Variant allele frequency")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 25)) 


ggsave(file="Figure_S7A.pdf",
       width=24,
       height=20,
       units="cm")


#-- Number of private and public targets, respectively
table(master.maf$label, useNA = "always")

#-- Number of public and private mutations
hehe<-unique(master.maf[,c("tumor", "Mutation", "label")])
table(hehe$label, useNA = "always")

#-- Private and Public mutations
master.sum <- master.sum %>% arrange(desc(mut.reads))

#-- Fraction of public mutations
master.sum <- master.sum %>% mutate(frac.public=mut.public/(mut.total-mut.unk))


#--- Compare public and private mutation MAFs
ind1<-which(master.maf$label=="private")
ind2<-which(master.maf$label=="public")

median(master.maf$MAF[ind1])
median(master.maf$MAF[ind2])

wilcox.test(master.maf$MAF[ind1], master.maf$MAF[ind2])

#-- Tail probabilities

print(sum(master.maf$MAF[ind1]<.4)/length(ind1))
print(sum(master.maf$MAF[ind2]<.4)/length(ind2))


#--------------------------------------------------
#--------------------------------------------------
#-- BRING IN DRIVER MUTATION ANNOTATION
#--------------------------------------------------
#--------------------------------------------------

#-- Convert master.maf to something more informative
master.maf$CHROM<-sub(".*chr", "", sub(":.*", "", master.maf$Mutation))
master.maf$POS<-sub("(.*):.*", "\\1", master.maf$Mutation)
master.maf$POS<-sub(".*:(.*)", "\\1", master.maf$POS)


MUT.annotation <- read.csv(file="Data/Annotations_SNV/SNV_annotations_complete.csv") %>%
  mutate(Driver=replace(Driver, Driver=="", "intergenic")) %>%
  rename(tumor=DCISsample)

#-- Render the tumor names compatible with each other
tumor_rename<- function(old) {
  new<-paste("Tumor_", sub(".*d", "", old), ".", sep="") 
  return(new)
}

MUT.annotation$tumor<-sapply(MUT.annotation$tumor, tumor_rename)

table(MUT.annotation$Driver)


#-- Update the master.maf with the original driver classification 

master.maf <- master.maf %>% 
  mutate(Driver =NA) %>%
  mutate(Gene=NA) %>%
  mutate(Ref=NA) %>%
  mutate(Alt=NA) %>%
  mutate(VariantType=NA)

for(k in 1: dim(master.maf)[1]){
  
  #-- Find the corresponding mutation for any given tumor/chr/position triple
  ind<-which(MUT.annotation$CHROM==master.maf$CHROM[k] &
               MUT.annotation$POS==master.maf$POS[k] &
               MUT.annotation$tumor==master.maf$tumor[k])
  
  #-- this would indicate an issue (can't find the mutation)
  if(length(ind)==0){print(k)}
  
  #-- update the driver status
  if(length(ind)>0){
    master.maf$Driver[k]<-MUT.annotation$Driver[ind]
    master.maf$Gene[k]<-MUT.annotation$GeneName[ind]
    master.maf$Ref[k]<-MUT.annotation$REF[ind]
    master.maf$Alt[k]<-MUT.annotation$ALT[ind]
    master.maf$VariantType[k]<-MUT.annotation$VariantType[ind]
  }
  
}

#-- Create a binary driver/passenger classification
master.maf <- master.maf %>%
  mutate(Driver.bin = ifelse(Driver=="likely", "Driver", "Passenger"))

#-- Get some useful statistics on the mutation level

he <- master.maf %>% select(-MAF)
he <- unique(he)

table(he$Driver.bin, useNA = "always")
table(he$label, useNA="always")


#-- Create a data frame with the driver mutations
driver.mutations <- master.maf %>% 
  filter(Driver.bin=="Driver") %>%
  select(-MAF) %>% distinct()

#-- Update master.sum sheet accordingly
master.sum <- master.sum %>% mutate(Drivers=0)

for(k in 1:dim(master.sum)[1]){
  ind<-which(driver.mutations$tumor ==master.sum$tumor[k])
  if(length(ind)>0){
    master.sum$Drivers[k]<-length(ind)
  }
}


#-------------------------------------
#-- Plot by detailed driver status
#-------------------------------------
ggplot(master.maf, aes(x=MAF, color=Driver, fill=Driver)) +
  geom_histogram(alpha=0.5, position="identity", binwidth = .03) +
  coord_cartesian(xlim=c(0,1)) + ggtitle("All Mutations and Spots") +
  facet_grid(~Driver)+theme_bw()


#-------------------------------------
#-- Plot by binary driver status (Figure S7B)
#-------------------------------------
ggplot(master.maf, aes(x=MAF,  fill=Driver.bin)) +
  scale_fill_manual(values=c("Driver"=pal1[4],"Passenger"=pal1[3] ))+
  geom_histogram(alpha=0.6, position="identity", binwidth = .025, mapping = aes(y = stat(ncount)))+
  coord_cartesian(xlim=c(0,1)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(text = element_text(size = 18)) +
  ylab("Count (normalized")+
  xlab("Minor allele frequency (MAF)")+
  theme(legend.position = "none")+
  theme(text = element_text(size = 25)) 


ggsave(file="Figure_S7B.pdf",
       width=24,
       height=20,
       units="cm")


dat.redox<-master.maf %>%
  select(tumor, CHROM, POS, Driver.bin, label)

dat.redox<-unique(dat.redox)

table(dat.redox$Driver.bin, dat.redox$label, useNA="always")

#-- Compare passenger and driver mutation MAFs
ind1<-which(master.maf$Driver.bin=="Driver")
ind2<-which(master.maf$Driver.bin=="Passenger")

median(master.maf$MAF[ind1])
median(master.maf$MAF[ind2])

wilcox.test(master.maf$MAF[ind1], master.maf$MAF[ind2])
ks.test(master.maf$MAF[ind1],master.maf$MAF[ind2])


#-------------------------------------
#-- Outputs
#-------------------------------------

nam.ord<-c("Tumor_03.",
           "Tumor_17.",
           "Tumor_24.",
           "Tumor_211.",
           "Tumor_213.",
           "Tumor_222.",
           "Tumor_225.",
           "Tumor_227.",
           "Tumor_315.",
           "Tumor_53.",
           "Tumor_65.",
           "Tumor_66.",
           "Tumor_88.",
           "Tumor_91.",
           "Tumor_168.",
           "Tumor_173.",
           "Tumor_191.",
           "Tumor_286.")

master.sum <- master.sum[match(nam.ord,master.sum$tumor),]
colSums(master.sum[,-1])


driver.mutations <- driver.mutations %>%
  select(tumor, CHROM, POS, Ref, Alt, Gene, VariantType) %>% 
  group_by(tumor)

#-- Update the tumor name format
driver.mutations$tumor<-sub("[.].*", "", driver.mutations$tumor)
driver.mutations$tumor<-as.numeric(as.character(sub(".*_", "", driver.mutations$tumor)))

#-- Update the tumor name format
master.sum$tumor<-sub("[.].*", "", master.sum$tumor)
master.sum$tumor<-as.numeric(as.character(sub(".*_", "", master.sum$tumor)))


#-- Mutations summary write the output
write.csv(master.sum, file=paste("Mutation_summary_", type, ".csv", sep=""), row.names = FALSE)


#-- Driver mutations write the output
write.csv(driver.mutations, file=paste("Driver_mutations_", type, ".csv", sep=""), row.names = FALSE)


#-- Summary for simulations data_M1_M2.csv
master.sum <- master.sum %>%
  select(-Drivers)

write.csv(master.sum, file="data_M1_M2.csv", row.names = FALSE)




