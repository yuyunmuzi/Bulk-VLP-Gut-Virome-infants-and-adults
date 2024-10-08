library(vegan)
library(ggplot2)
library(vegan)
library(ggpubr) 
library(scales) 
library(ggstatsplot)
library(dplyr) 
library(reshape)
library(tidyr) 
library(hrbrthemes)
library(viridis) 
library(ggalluvial)
library(hrbrthemes)
library(stringr)
library(rlang)
library(ggbreak)
library(patchwork)


#PRJNA524703,infant143--------
###PRJNA835720 for adults (151) is similar to infants ----
load("./AdditionalData/infant/infant.bulk_1.RData")
load("./AdditionalData/infant/infant.vlp_1.RData")
load("./AdditionalData/adult/adult.bulk_1.RData")
load("./AdditionalData/adult/adult.vlp_1.RData")
#1. shannon index ------------
R1 <- intersect(rownames(NC.bac11),rownames(NC.mNGS_1))
NC.bac2 <- subset(NC.bac11,rownames(NC.bac11) %in% R1)
NC.mNGS_2 <- subset(NC.mNGS_1,rownames(NC.mNGS_1) %in% R1)

b_shannon <- diversity(NC.bac2,index="shannon")
m_shannon <- diversity(NC.mNGS_2,index="shannon")
correlation<-cbind(b_shannon,v_shannon) %>% as.data.frame()
cor.test(b_shannon,v_shannon, alternative = "two.sided", method = "spearman",conf.level = 0.95)
infantp1_shannon <- ggplot(correlation,aes(x=b_shannon,y=m_shannon))+
    stat_smooth(method=lm,formula=y~x,color="black")+
    geom_point(color="#CB7368")+
    theme_bw()+
    annotate('text', label = 'bulk\nrho = ##\np-value = ## ',
             x = 0, y =4, size = 4,hjust = 0)+
    xlab("bacteriome shannon")+
    ylab("virome shannon")

#richness----------
b_ob<-rowSums( NC.bac2!=0)
v_ob<-rowSums(NC.mNGS_2!=0)
correlation<-cbind(b_ob,v_ob)
correlation<-as.data.frame(correlation)
library(ggplot2)
ggplot(correlation,aes(x=b_ob,y=v_ob))+
  stat_smooth(method=lm,formula=y~x,color="black")+
  geom_point(color="#CB7368")+
  theme_bw()+
  xlab("bacteriome richness")+
  ylab("virome richness")+
  annotate('text', label = 'Bulk\nrho = ##\np-value = ## ',
           x = 2, y =220, size = 4,hjust = 0) 

#2.overlap------------------
NC.m.overlap <- read.delim(file="F:/03噬菌体/2024NC_infant/NC.m_speciesp.contig",header =FALSE, sep = "\t")
NC.mNGS_3 <- NC.mNGS_2[,colnames(NC.mNGS_2) %in% NC.m.overlap$V1]
NC.mNGS_3 <- na.omit(NC.mNGS_3)

b_shannon <- diversity(NC.bac2,index="shannon")
m_shannon <- diversity(NC.mNGS_3,index="shannon")
correlation<-cbind(b_shannon,m_shannon)
correlation<-as.data.frame(correlation)
cor.test(b_shannon,m_shannon, alternative = "two.sided", method = "spearman",conf.level = 0.95)
  infantp1_shannon <- ggplot(correlation,aes(x=b_shannon,y=m_shannon))+
    stat_smooth(method=lm,formula=y~x,color="black")+
    geom_point(color="#CB7368")+
    theme_bw()+
    annotate('text', label = 'bulk\nrho = ##\np-value = ## ',
             x = 0, y =4, size = 4,hjust = 0)+
    xlab("bacteriome shannon")+
    ylab("virome shannon")
  
b_ob<-rowSums( NC.bac2!=0)
v_ob<-rowSums(NC.mNGS_2!=0)
correlation<-cbind(b_ob,v_ob)
correlation<-as.data.frame(correlation)
cor.test(b_ob,v_ob, alternative = "two.sided", method = "spearman",conf.level = 0.95)
infantp1_overm <-  ggplot(correlation,aes(x=b_ob,y=v_ob))+
  stat_smooth(method=lm,formula=y~x,color="black")+
  geom_point(color="#CB7368")+
  theme_bw()+
  xlab("bacteriome richness")+
  ylab("virome richness")+
  annotate('text', label = 'Bulk\nrho = ##\np-value = ## ',
           x = 2, y =140, size = 4,hjust = 0) 



#3.Accumulation curve--------------
rr <- intersect(rownames(NC.mNGS_1),rownames(NC.vlp_1))
NC.mNGS_1.count <- subset(NC.mNGS_1,rownames(NC.mNGS_1) %in% rr)
NC.vlp_1.count <- subset(NC.vlp_1,rownames(NC.vlp_1) %in% rr)

contig_unique="SRR8652857_k141_1055"
contig_count<-data.frame(matrix(0, nrow = 150, ncol = 2))

for (iii in 1:141) {
  a<-colnames(NC.mNGS_1.count[,which(NC.mNGS_1.count[iii,]!=0)])
  contig_unique=unique(c(contig_unique,unique(a)))
  contig_count[iii,1]=length(contig_unique)
  contig_count[iii,2]=iii
}
contig_unique="SRR8652812_k119_111"
contig_count1<-data.frame(matrix(0, nrow = 138, ncol = 2))
for (iii in 1:152) {
  a<-colnames(NC.vlp_1.count[,which(NC.vlp_1.count[iii,]!=0)])
  contig_unique=unique(c(contig_unique,unique(a)))
  contig_count1[iii,1]=length(contig_unique)
  contig_count1[iii,2]=iii
}
"mNGS"="#CB7368","vNGS"="#75AADB"
NC.p1 <- ggplot() + 
  geom_line(data = contig_count, aes(x = X2, y = X1), color = "#CB7368", size = 1.5) +
  geom_line(data = contig_count1, aes(x = X2, y = X1), color = "#75AADB", size = 1.5) +
  xlab("Sample") +
  ylab("contig") +
  theme_bw()+
  scale_y_continuous(breaks =seq( 0, max(contig_count1$X1), by = 2000))
NC.p1+ad.p1

#4. abundance----------------

##碱基----
infant_v_bases<- read.delim(file="F:/03噬菌体/2024NC_infant/infant_v_bases.csv",header =FALSE, sep = "\t")
infant_v_bases$V2 <- as.numeric(gsub(",", "", infant_v_bases$V2))
v_observe<-data.frame(rowSums(NC.vlp_1!=0))
colnames(v_observe)="richness"
v_observe$group <- rownames(v_observe)
v_observe <- merge(v_observe,infant_v_bases,by.x="group",by.y="V1")
v_observe <- na.omit(v_observe)
v_observe$percent <- v_observe$richness/v_observe$V2*1000000

infant_m_bases<- read.delim(file="F:/03噬菌体/2024NC_infant/infant_m_bases.csv",header =FALSE, sep = "\t")
infant_m_bases$V2 <- as.numeric(gsub(",", "", infant_m_bases$V2))
infant_m_bases <- merge(infant_m_bases,m.tov,by.x="V1",by.y="group")
infant_m_bases <- infant_m_bases[,c(2,3)]
m_observe<-data.frame(rowSums(NC.mNGS_1!=0))
m_observe <- na.omit(m_observe)
colnames(m_observe)="richness"
m_observe$group <- rownames(m_observe)
m_observe <- merge(m_observe,infant_m_bases,by.x="group",by.y="v_group")
m_observe <- na.omit(m_observe)
m_observe$percent <- m_observe$richness/m_observe$V2*1000000

a<-data.frame(group="mNGS",rich_base=m_observe$V2)
d<-data.frame(group="vNGS",rich_base=v_observe$V2)
all<-rbind(a,d)
NC.base1 <- all %>% 
  mutate(group=factor(group,levels=c("mNGS","vNGS"))) %>% 
  ggplot(aes(group,rich_base/1000000,color=group)) +
  geom_boxplot(width=0.3,size=1) +scale_y_log10()+
  scale_color_manual(values =c("mNGS"="#CB7368","vNGS"="#75AADB")) +
  #geom_jitter(shape=16, position = position_jitter(0.2))+
  theme_bw()+ #换个主题
  ylab("Nr. of million reads")+
  stat_compare_means(comparisons = list(c("mNGS","vNGS") ),
                     label = "p.signif",
                     method = "wilcox.test",
                     paired = TRUE)

###Relative abundance of vOTUs----
median.na <- function(x){y <- median(x[x>0],na.rm = T); return(y)}
a<-data.frame(group="mNGS",abun=apply(NC.mNGS_3,2,median.na))
d<-data.frame(group="vNGS",abun=apply(NC.vlp_3,2,median.na))
all<-rbind(a,d)
all$abun <- all$abun*100
all <- all[which(all$abun<5),]
all %>% ggplot(aes(x = abun, fill = group)) +
  geom_density(alpha = 0.7) +
  scale_fill_manual(values = c("mNGS"="#CB7368","vNGS"="#75AADB")) +
  theme_bw() +
  scale_x_log10()+
  geom_vline(data=all %>% group_by(group) %>% summarise(median = median((abun))),
             aes(xintercept = median, color = group), linetype = "dashed",size=0.75) +
  labs(x = "Abundance", y = "Density") +
  scale_color_manual(values = c("mNGS"="#CB7368","vNGS"="#75AADB"))


##The dominant vOTUs-----
m_abun <- data.frame(abun=c(1:143))
for(i in 1:143){
  m_abun$abun[i] <- max(NC.mNGS_1[i,])
}
m_abun$group="mNGS"
m_abun$Group <- rownames(NC.mNGS_1)
col_names <- rownames(NC.mNGS_1)

v_abun <- data.frame(abun=c(1:141))
for(i in 1:141){
  v_abun$abun[i] <- max(NC.vlp_1[i,])
}
v_abun$group="vNGS"
v_abun$Group <- rownames(NC.vlp_1)

col_names <- colnames(NC.vlp_1)
max_col_names <- col_names[apply(NC.vlp_1, 1, which.max)]
v_abun$col_name <- max_col_names

high_abun <- rbind(m_abun,v_abun)
high_abun1 <- merge(high_abun,NC.month,by.x="Group",by.y="group")
high_abun1 <- na.omit(high_abun1)
high_abun1 %>% 
  ggplot(aes(group,abun,color=group)) +
  geom_boxplot(width=0.3,size=1) +
  scale_color_manual(values =c("mNGS"="#CB7368","vNGS"="#75AADB")) +
  theme_bw()+  
  ylab("Highest abundance in each sample")+
  stat_compare_means(comparisons = compaired ,
                     label = "p.signif",
                     method = "wilcox.test",
                     paired = TRUE)+
  facet_grid(~ factor(month, levels = c("Month 0","Month 1","Month 4","Adult")))

#5. Comparison of the alpha diversity of the bacteriome and virome in the samples
####richness----
NC.month<-read.delim(file="F:/03噬菌体/2024NC_infant/group_month.txt",header =TRUE, sep = "\t")
NC.m_ob<-data.frame(rowSums(NC.mNGS_1!=0))
colnames(NC.m_ob)="richnness"
NC.m_ob$group=rownames(NC.m_ob)
NC.m_ob <- na.omit(NC.m_ob)
NC.v_ob<-data.frame(rowSums(NC.vlp_1!=0))
colnames(NC.v_ob)="richnness"
NC.v_ob$group=rownames(NC.v_ob)
NC.v_ob <- na.omit(NC.v_ob)
NC.ob <- left_join(NC.m_ob,NC.v_ob,by="group")
NC.ob <- na.omit(NC.ob)
NC.ob <- left_join(NC.ob,NC.month,by="group")

colnames(NC.ob)=c("m_richness","group","v_richness","month")
NC.ob <- NC.ob[,c("m_richness","v_richness","group","month")]
NC.ob_1 <- data.frame(group=c(rep("mNGS",134),rep("vNGS",134)),observed_species=c(NC.ob$m_richness,NC.ob$v_richness),month=NC.ob$month)
NC.p2 <- ggplot(NC.ob_1,aes(group,observed_species,fill=group)) +
  geom_boxplot(width=0.3,size=1) +
  scale_fill_manual(values =c("mNGS"="#CB7368","vNGS"="#75AADB")) +
  theme_bw()+  
  stat_compare_means(comparisons = list(c("mNGS","vNGS")),
                     label = "p.signif",
                     method = "wilcox.test",
                     paired = TRUE)+
  facet_grid(~ factor(month, levels = c("Month 0","Month 1","Month 4")))

#shannon----
mShannon_phage <- diversity(NC.mNGS_1,index="shannon")
VShannon_phage <- diversity(NC.vlp_1,index="shannon")
vNGS<-data.frame(group=rep("mNGS",143),shannon=mShannon_phage)
vNGS<-data.frame(group=rep("vNGS",141),shannon=VShannon_phage)
mmvv<-rbind(mNGS,vNGS)
library(ggpubr)
compaired <- list(c("mNGS","vNGS"))
shannon_infant <- ggplot(mmvv,aes(group,shannon,color=group)) +
  geom_boxplot(width=0.3,size=1) +
  scale_color_manual(values =c("mNGS"="#75AADB","vNGS"="#E89138")) +
  theme_bw()+  
  stat_compare_means(comparisons = compaired ,
                     label = "p.signif",
                     method = "wilcox.test",
                     paired = TRUE)

####evenness----
mmvv<-cbind(NC.shannon_1,NC.ob_1[,2])
mmvv$eveness <- mmvv$shannon/log2(mmvv$`NC.ob_1[, 2]`)
compaired <- list(c("mNGS","vNGS"))
even.pp1 <- mmvv %>% 
  ggplot(aes(group,eveness,color=group)) +
  geom_boxplot(width=0.3,size=1) +
  scale_color_manual(values =c("mNGS"="#CB7368","vNGS"="#75AADB")) +
  theme_bw()+  
  stat_compare_means(comparisons = compaired ,
                                label = "p.signif",
                     method = "wilcox.test",paired = TRUE)

####beta-diveisty----
NC.mNGS_1_0 <- subset(NC.mNGS_1,rownames(NC.mNGS_1) %in% NC.month[NC.month$month=="Month 0",1])
NC.mNGS_1_1 <- subset(NC.mNGS_1,rownames(NC.mNGS_1) %in% NC.month[NC.month$month=="Month 1",1])
NC.mNGS_1_4 <- subset(NC.mNGS_1,rownames(NC.mNGS_1) %in% NC.month[NC.month$month=="Month 4",1])

NC.m_bd0<-vegdist(na.omit(NC.mNGS_1_0),method = "bray")
NC.m_bd1<-vegdist(na.omit(NC.mNGS_1_1),method = "bray")
NC.m_bd4<-vegdist(na.omit(NC.mNGS_1_4),method = "bray")

NC.vNGS_1_0 <- subset(NC.vlp_1,rownames(NC.vlp_1) %in% NC.month[NC.month$month=="Month 0",1])
NC.vNGS_1_1 <- subset(NC.vlp_1,rownames(NC.vlp_1) %in% NC.month[NC.month$month=="Month 1",1])
NC.vNGS_1_4 <- subset(NC.vlp_1,rownames(NC.vlp_1) %in% NC.month[NC.month$month=="Month 4",1])

NC.v_bd0<-vegdist(na.omit(NC.vNGS_1_0),method = "bray")
NC.v_bd1<-vegdist(na.omit(NC.vNGS_1_1),method = "bray")
NC.v_bd4<-vegdist(na.omit(NC.vNGS_1_4),method = "bray")

bary <- data.frame(group=c(rep("mNGS",5245),rep("vNGS",5596)),bray_cruit=c(NC.m_bd0,NC.m_bd1,NC.m_bd4,NC.v_bd0,NC.v_bd1,NC.v_bd4),
                   month=c(rep("Month 0",105),rep("Month 1",190),rep("Month 4",4950),
                           rep("Month 0",153),rep("Month 1",190),rep("Month 4",5253)))
bary <- bary[which(bary$bray_cruit>0.5),]
NC.pp1 <- ggplot(bary,aes(group,bray_cruit,fill=group)) +
  geom_boxplot(width=0.3,size=0.75) +
  scale_fill_manual(values =c("mNGS"="#CB7368","vNGS"="#75AADB")) +
  theme_bw()+ 
  stat_compare_means(comparisons = list(c("mNGS","vNGS")),
                     label = "p.signif",
                     method = "wilcox.test",paired = TRUE)+
  facet_grid(~ factor(month, levels = c("Month 0","Month 1","Month 4")))

##co-identity---
load("./AdditionalData/infant/infant.bulk_1a1.RData")
load("./AdditionalData/infant/infant.vlp_1a1.RData")
load("./AdditionalData/adult/adult.bulk_1a1.RData")
load("./AdditionalData/adult/adult.vlp_1a1.RData")
NC.vlp_1a2 <- NC.vlp_1a1 
  NC.vlp_1a2[NC.vlp_1a2>0]=1
  NC.mNGS_1a2 <- NC.mNGS_1a1 
  NC.mNGS_1a2[NC.mNGS_1a2>0]=1
  R <- intersect(colnames(NC.vlp_1a2),colnames(NC.mNGS_1a2))
  aa <- m_prevalence
  aa$contribution <- numeric(nrow(m_prevalence))
  for(i in 1:length(m_prevalence$Genus)){
    a=0
    if (m_prevalence$Genus[i] %in% rownames(NC.vlp_1a2)) {
      for(j in R){
        if(NC.mNGS_1a2[rownames(NC.mNGS_1a2) == m_prevalence$Genus[i], j]==1 && NC.vlp_1a2[rownames(NC.vlp_1a2) == m_prevalence$Genus[i], j]==1 ){
          a=a+1
        }
      }
      m=NC.mNGS_1a2[rownames(NC.mNGS_1a2) == m_prevalence$Genus[i], ]
      m <- names(m)[m > 0]
      v=NC.vlp_1a2[rownames(NC.vlp_1a2) == m_prevalence$Genus[i], ]
      v <- names(v)[v > 0]
      l=length(unique(m,v))
    }
      aa$contribution[i]=round(a/l,3)
  }
  
  library(ggplot2)
  library(dplyr)
  
  # The dataset is provided in the gapminder library
  #library(gapminder)
  #data <- gapminder %>% filter(year=="2007") %>% dplyr::select(-year)
  
  aa$group="mNGS"
  bb <- left_join(v_prevalence,aa,by="Genus")
  bb <- na.omit(bb)
  #bb <- bb[which(bb$contribution>0),]
  bb <- bb[c(1,2,4)]
  bb$group="vNGS"
  colnames(bb)[1]="prevalence"
  #aa1 <- aa[which(aa$contribution>0),]
  A <- rbind(aa,bb)
  
  A %>%
    arrange(desc(prevalence)) %>%
   # mutate(country = factor(country, country)) %>%
    ggplot(aes(x=prevalence, y=reorder(Genus,prevalence), size = prevalence,color=contribution)) +
    geom_point(alpha=0.5) +
    scale_color_gradient(low = "lightblue", high = "#142D45")+
    facet_wrap(~ group,scales = "free_y",ncol = 2)+ # scales = "free_y", 
    #scale_size(range = c(.1, 20), name="Prevalence")+
    theme_classic()+
    ylab("Genus (present in both mNGS and vNGS)")+
    xlab("Prevalence") +
    scale_x_continuous(limits = c(0, 1),labels = percent_format()) 
  cor.test(aa$prevalence,aa$contribution, alternative = "two.sided", method = "spearman",conf.level = 0.95)
  p1 <- aa %>%  ggplot(aes(x=prevalence,y=contribution))+
   stat_smooth(method="loess",formula=y~x,color="#CB7368")+
   geom_point(color="#CB7368")+
   theme_bw()+ylim(c(0,1))+
    annotate('text', label = 'rho = 0.44\np-value = 1.335e-15-05 ',
             x = 0.5, y =0.9, size = 4,hjust = 0)+
    xlab("Prevalence of genus")+
    ylab("The consistent of genus of bulk and VLP")
  cor.test(bb$prevalence,bb$contribution, alternative = "two.sided", method = "spearman",conf.level = 0.95)
 p2 <-  bb %>%  ggplot(aes(x=prevalence,y=contribution))+
    stat_smooth(method="loess",formula=y~x,color="#75AADB")+
    geom_point(color="#75AADB")+
    theme_bw()+ylim(c(0,1))+
    annotate('text', label = 'rho = 0.53\np-value < 2.2e-16 ',
             x = 0.4, y =0.9, size = 4,hjust = 0)+
    xlab("Prevalence of genus")+
    ylab("The consistent of genus of bulk and VLP")
 p1+p2 

##LEfse----
colnames(our.mNGS_1a2) <- paste(colnames(our.mNGS_1a2), "_mNGS", sep = "")
colnames(our.vlp_1a1) <- paste(colnames(our.vlp_1a1), "_vlp", sep = "")

our.mNGS_1a2 <- our.mNGS_1a2[!rownames(our.mNGS_1a2) %in% E_genus.list$V1, ]#433 150
our.vlp_1a1 <- our.vlp_1a1[!rownames(our.vlp_1a1) %in% E_genus.list$V1, ]#459  150

our.mNGS_1a2 <- data.frame(our.mNGS_1a2[,apply(our.mNGS_1a2,2,max)>0]) #190 150
our.vlp_1a1 <- data.frame(our.vlp_1a1[,apply(our.vlp_1a1,2,max)>0])  #183  150
feat_list <- list()
feat_list[["our.mNGS_1a2"]] <- as.data.frame(t(our.mNGS_1a2))  #150  190
feat_list[["our.vlp_1a1"]] <- as.data.frame(t(our.vlp_1a1))   #150  183
my_pair.table=function(feat_list){
  library(tidyverse)
  
  
  feat_name=lapply(feat_list,function(data){colnames(data)})
  feat_ID=feat_name%>%unlist()%>%as.vector()%>%unique()
  b=lapply(feat_list, function(data){
    add=setdiff(feat_ID,colnames(data))
  })
  
  my_add=function(b,feat){
    if(length(b)!=0){
      b=as.character(b)
      x_add=data.frame(matrix(0,dim(feat)[1],length(b)))
      colnames(x_add)=b
      xtest=cbind(feat,x_add)
    }else{
      xtest=feat
    }
    return(xtest)
  }
  
  data_add=list()
  for (i in 1:length(feat_list)) {
    data_add[[i]]=my_add(b[[i]],feat_list[[i]])
  }
  names(data_add)=names(feat_list)
  
  feat_seq=sort(feat_ID)
  data_add=lapply(data_add, function(data){new_data=data[,c(feat_seq)]})
  
  return(data_add)
}
feat_list <- my_pair.table(feat_list) #150  190/150  183

our <- rbind(feat_list$our.mNGS_1a2,feat_list$our.vlp_1a1)

#
feat_list <- list()
feat_list[["our"]] <- as.data.frame(our)#300 261
rownames_meta <- rownames(our)
group_meta <- ifelse(grepl("_mNGS$", rownames_meta), "Control", "Case")
meta_our <- data.frame(sample_id = rownames_meta, Group = group_meta,Bodysite="oral",disease_stage=NA,country="us",sex=NA,host_age=NA,BMI=NA)
meta_our <- column_to_rownames(meta_our,var = "sample_id")

meta_list <- list()
meta_list[["our"]] <- meta_our

source("./AdditionalData/siamcat_models_adj.R")
# library(pkgbuild)
# library(profvis)
# library(devtools)
# library(MMUPHin)
# devtools::install_github("rstudio/r2d3")
caries_mymarker.adj <- my_marker_adj(feat_list,meta_list,NULL,NULL,T,T,T,
                                     is_plot=F,lda_cutoff=3,nproj_cutoff=1,level="genus",change_name=F)


caries_marker <- caries_mymarker.adj$marker_data
caries_marker$scientific_name <- gsub("g__","", caries_marker$scientific_name)
caries_marker$LDA <- as.numeric(caries_marker$LDA)
caries_marker$p_value <- as.numeric(caries_marker$p_value)
summary(caries_marker)
sum(caries_marker$class == "adjust")
sum(caries_marker$p_value <= 0.05)

caries_marker_filter <- caries_marker %>%
  filter(abs(LDA) >= 3,p_value <= 0.05) %>% 
  filter(class == "adjust")
unique(caries_marker_filter$scientific_name) #146个

caries_marker_filter <- caries_marker_filter %>%
  mutate(LDA_sign = ifelse(LDA > 0, 'vNGS', 'mNGS'))

LDA_sign_changes <- caries_marker_filter %>%
  group_by(scientific_name) %>%
  summarise(Changes = n_distinct(LDA_sign)) %>%
  filter(Changes == 1) 

caries_marker_filter_stable <- caries_marker_filter %>%
  semi_join(LDA_sign_changes, by = 'scientific_name')
unique(caries_marker_filter_stable$scientific_name) #139


negative_names <- caries_marker_filter_stable %>%
  filter(LDA_sign == "mNGS") %>%
  arrange(LDA) %>%
  pull(scientific_name)

positive_names <- caries_marker_filter_stable %>%
  filter(LDA_sign == "vNGS") %>%
  arrange(LDA) %>%
  pull(scientific_name)

ordered_names <- c(negative_names, positive_names)

caries_marker_filter_stable$scientific_name <- 
  factor(caries_marker_filter_stable$scientific_name, levels = unique(ordered_names))
ggplot(caries_marker_filter_stable, aes(x = LDA, y = scientific_name,fill= LDA_sign)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#C45C69","#4D779B"),name = "LDA_sign")+
  scale_x_continuous(breaks = c(seq(-4, -2, by = 1), seq(2, 4, by = 1))) +
  #scale_x_continuous(breaks = c(seq(min(caries_marker_filter_stable$LDA), -2, by = 1), seq(2, max(caries_marker_filter_stable$LDA), by = 1))) +
  theme_classic() +
  labs(x = "LDA", y = "Scientific Name") +
  theme(axis.text.y = element_text(angle = 45, hjust = 1))
