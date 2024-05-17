library(tidyverse)
library(reshape)
library(tidyr) 

#PRJNA524703,infant141--------
###PRJNA835720 for adults (151) is similar to infants ----

NC.vlp<-read.delim(file="F:/2024NC_infant/NC.all.rpkm75.txt",header =FALSE, sep = "\t")
NC.mNGS<-read.delim(file="F:/2024NC_infant/mNGS.infant.coverm.txt",header =FALSE, sep = "\t")
NC.mNGS2<-read.delim(file="F:/2024NC_infant/01_mNGS.rpkm75.infant.txt",header =FALSE, sep = "\t")
NC.bac <- read.delim(file="F:/2024NC_infant/merged_abundance_species.txt",header =TRUE, sep = "\t")

####PRJNA524703 vlp group bac----
NC.bac1 <- data.frame(t(NC.bac))
colnames(NC.bac1)=NC.bac1[1,]
NC.bac1 <- NC.bac1[-1,]
NC.bac1$group <- rownames(NC.bac1)
NC.bac1 <- left_join(NC.bac1,m.tov,by="group")
rownames(NC.bac1)=NC.bac1$v_group
NC.bac1 <- NC.bac1[,c(-533,-532)]
bac.names=rownames(NC.bac1)
NC.bac1<-as.data.frame(lapply(NC.bac1,as.numeric))
NC.bac1<-sweep(NC.bac1,1,rowSums(NC.bac1),`/`)
rownames(NC.bac1)=bac.names
m.tov <- read.delim(file="F:/2024NC_infant/m_to_v.txt",header =FALSE, sep = "\t")
colnames(m.tov)=c("group","v_group")
NC.bac11 <- NC.bac1[rowSums(is.na(NC.bac1)) != ncol(NC.bac1), ]

#NC.mNGS_1,PRJNA524703,infant141mNGS 95not remove----------
colnames(NC.mNGS)=c("group","contig","Abuandance")
NC.mNGS_1<-spread(NC.mNGS, group, Abuandance)
NC.mNGS_1[is.na(NC.mNGS_1)] <- 0
temp=NC.mNGS_1[,1]
NC.mNGS_1=NC.mNGS_1[,-1]
rownames(NC.mNGS_1)=temp
NC.rm.mNGS <- colnames(NC.mNGS_1)
NC.mNGS_1<-as.data.frame(t(NC.mNGS_1))
NC.mNGS_1<-as.data.frame(lapply(NC.mNGS_1,as.numeric))
NC.mNGS_1<-sweep(NC.mNGS_1,1,rowSums(NC.mNGS_1),`/`)
rownames(NC.mNGS_1) <- NC.rm.mNGS
NC.mNGS_1$group <- rownames(NC.mNGS_1)
NC.mNGS_1 <- left_join(NC.mNGS_1,m.tov,by="group")
rownames(NC.mNGS_1)=NC.mNGS_1$v_group
NC.mNGS_1 <- NC.mNGS_1[,c(-4110,-4109)]

####NC.vlp_1,PRJNA524703 infant 141VLP not 95remove-----
NC.vlp<-NC.vlp[,-3]
colnames(NC.vlp)=c("group","contig","Abuandance")
rm(NC.vlp_1)
NC.vlp <- unique(NC.vlp)
NC.vlp_1<-spread(NC.vlp, group, Abuandance)
NC.vlp_1[is.na(NC.vlp_1)] <- 0
temp=NC.vlp_1[,1]
NC.vlp_1=NC.vlp_1[,-1]
rownames(NC.vlp_1)=temp
NC.rm <- colnames(NC.vlp_1)
NC.vlp_1<-as.data.frame(t(NC.vlp_1))
NC.vlp_1<-as.data.frame(lapply(NC.vlp_1,as.numeric))
NC.vlp_1<-sweep(NC.vlp_1,1,rowSums(NC.vlp_1),`/`)
rownames(NC.vlp_1) <- NC.rm



