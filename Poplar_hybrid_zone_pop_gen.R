
#### Bolte et al., 2024
#### Populus trichocarpa x balsamifera hybrid zones
## Dr. Connie Bolte, conniebolte@gmail.com

#################### PCA #####################

###PCA for 012 coded  vcf files
##Following method in Patterson et al 2006

library(data.table)
install.packages("tidyverse")
library(tidyverse)
install.packages("ggsci")
library(ggsci)
install.packages("ggforce")
library(ggforce)
install.packages("devtools")
library(devtools)
install.packages("viridis")
library(viridis)
#install.packages("LEA")
#library(LEA)
install.packages("lme4")
library(lme4)
install.packages("lmerTest")
library(lmerTest)
library(readr)
install.packages("plyr")
library(plyr)
library(readxl)
install.packages("wesanderson")
library(wesanderson)

#source('~/Desktop/Imports.R')

##################### PCA_gen ##########################

##### input files #####

### df_gen: genotypic data with individuals as rows and snps as coloumns.
###         Can include missing data. Either genotype probabilities or 012 format

### indv:  dataframe with rows corresponding to individuals in df_gen file 
###        Must have Pop and ID coloumn 

##### output files ####

### pca_out: 

### $ 'pca_df': dataframe with rows as individuals and coloumns as PC1-30, Pop, ID
### $ 'pve': list of proportion of variance explained for each PC


##### function #######

PCA_gen <- function(df_gen,indv,num=5){ #add ggplot, add tw, add # of output
  #pkgTest("ggplot2")
  
  df_gen <- apply(df_gen, 2, function(df) gsub(-1,NA,df,fixed=TRUE))
  df_gen <- apply(df_gen, 2, function(df) as.numeric(df))
  
  colmean<-apply(df_gen,2,mean,na.rm=TRUE)
  
  normalize<-matrix(nrow = nrow(df_gen),ncol=ncol(df_gen))
  af<-colmean/2
  
  for (m in 1:length(af)){
    nr<-df_gen[ ,m]-colmean[m]
    dn<-sqrt(af[m]*(1-af[m]))
    normalize[ ,m]<-nr/dn
  }
  
  normalize[is.na(normalize)]<-0
  
  method1<-prcomp(normalize, scale. = FALSE,center = FALSE)
  pve <- summary(method1)$importance[2,]
  print(pve[1:300])
  
  if(nrow(df_gen) < num){
    num <- nrow(df_gen)
  }
  
  pca_X<- method1$x[,1:num]
  
  pca_X <- as.data.frame(pca_X)
  pca_X$Pop <- indv$Pop
  pca_X$ID <- indv$ID
  
  pca_out <- list("pca_df"=pca_X,"pve"=pve)
  
  #print(PCA_fig(pca_out))
  
  return(pca_out)
  return(method1)
}





##### PCA of merged angustifolia samples with the 575 sample collected for this hybrid zone #####

all_samps <- read.table("~/Documents/Paper1/data_files/merged_angust_tri_bal_poplars_imputed_no_scaffolds_mac2_miss0.3.recode.vcf.012.indv")
head(all_samps)
all_samps_012 <- read.table("~/Documents/Paper1/data_files/merged_angust_tri_bal_poplars_imputed_no_scaffolds_mac2_miss0.3.recode.vcf.012")
View(all_samps_012[1:10,1:10])
df012 <- all_samps_012[,-1]


df012_na <- apply(df012, 2, function(df) gsub('9','NA',df,fixed=TRUE)) # no missing data so this function is obsolete.
View(df012_na[1:20,1:20])
formated<-apply(df012_na,2,function(df) as.numeric(df))
View(formated[1:20, 1:20])

colnames(all_samps) <- "ID"

K4_scores_623_pop_ids <- read.csv("~/Documents/Paper1/data_files/K4_scores_623_sample_pop_ids_angust_imputed_data.csv")
Pop <- K4_scores_623_pop_ids$Pop
pops_623 <- left_join(all_samps, K4_scores_623_pop_ids, by="ID")
head(pops_623)
pops_623 <- data.frame(ID=pops_623$ID, Pop=pops_623$Pop)


pca_out <- PCA_gen(formated,pops_623)

head(pca_out)
pve <- pca_out$pve[1:5]
pve
#    PC1     PC2     PC3     PC4     PC5 
# 0.22593 0.10657 0.05061 0.03896 0.02160 

pca_df <- pca_out$pca_df[,1:5]
pca_df <- cbind(pops_623,pca_df)

write.csv(pca_df,'~/Documents/Paper1/pca_df_623_samples_with_merged_angust_imputed_data.csv',row.names = FALSE)
pca_df <- read.csv("~/Documents/Paper1/pca_df_623_samples_with_merged_angust_imputed_data.csv")


############################ PLOT gradient color by Q-score ################################

pve <- c(0.22593, 0.10657, 0.05061, 0.03896, 0.02160 )

col= c("orange", "black", "green3")
col2= c("orange", "red", "lightblue")
pca_623 <- ggplot(data = pca_df, aes(x=PC1,y=PC2, fill=Pop)) + 
  geom_point(colour='black', size=3, pch=21) +
  scale_fill_manual(values=col2) +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))+
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=8), 
        axis.title = element_text(size = 10, colour="black",face = "bold",vjust = 1),
        #legend.text = element_text(size=16),
        #legend.title = element_text(size = 16, colour="black", face="bold", vjust=1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(pca_623)


ggsave('~/Documents/Paper1/PCA12_623_samples_with_merged_imputed_data_v2.pdf',pca_623,height=4,width=6,units='in')



pca_623b <- ggplot(data = pca_df, aes(x=PC2,y=PC3, fill=Pop)) + 
  geom_point(colour='black', size=3, pch=21) +
  scale_fill_manual(values=col2) +
  xlab(paste("PC",2," (",pve[2]*100,"%)",sep="")) + ylab(paste("PC",3," (",pve[3]*100,"%)",sep=""))+
  theme_bw() + 
  theme(legend.position = 'none',
        axis.text = element_text(size=8), 
        axis.title = element_text(size = 10, colour="black",face = "bold",vjust = 1),
        #legend.text = element_text(size=16),
        #legend.title = element_text(size = 16, colour="black", face="bold", vjust=1),
        panel.border = element_rect(size = 1.5, colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())


plot(pca_623b)


ggsave('~/Documents/Paper1/PCA23_623_samples_with_merged_imputed_data_v2.pdf',pca_623b,height=4,width=6,units='in')



pca_623c <- ggplot(data = pca_df, aes(x=PC3,y=PC4, fill= Pop)) + 
  geom_point(colour='black', size=3, pch=21) +
  scale_fill_manual(values=col2) +
  xlab(paste("PC",3," (",pve[3]*100,"%)",sep="")) + ylab(paste("PC",4," (",pve[4]*100,"%)",sep=""))+
  theme_bw() + 
  theme(#legend.position = 'none',
    axis.text = element_text(size=8), 
    axis.title = element_text(size = 10, colour="black",face = "bold",vjust = 1),
    legend.text = element_text(size=14),
    legend.title = element_text(size = 14, colour="black", face="bold", vjust=1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())


plot(pca_623c)


ggsave('~/Documents/Paper1/PCA34_623_samples_with_merged_imputed_data_v2.pdf',pca_623c,height=4,width=6,units='in')



full <- pca_623 + pca_623b + pca_623c + plot_layout(ncol=3)

ggsave('~/Documents/Paper1/PCA_triplot_PC12_PC23_pC34_623_samples_with_merged_imputed_data_v2.pdf',full,height=3,width=9,units='in')








##################### PCA of bi-allelic, LD pruned, SNP data for 546 samples. ##########################


df012<-fread("~/Documents/Paper1/data_files/poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012",sep="\t", data.table=F) 
df012 <- df012[,-1]
dim(df012)


df012_na <- apply(df012, 2, function(df) gsub('-1','NA',df,fixed=TRUE))
View(df012_na[1:20,1:20])
formated<-apply(df012_na,2,function(df) as.numeric(df))
View(formated[1:20, 1:20])
indv_LDprune01 <- read.table("~/Documents/Paper1/data_files/poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012.indv")
colnames(indv_LDprune01) <- "ID"


ID_clim_geo_transect <- read.csv("~/Documents/Paper1/data_files/Modified_transect_assignments_with_ClimGeo_data.csv")
ID_transect <- data.frame(ID=ID_clim_geo_transect$ID, Pop=ID_clim_geo_transect$Transect_SL)
Pop_ID_Sum <- ID_transect

pca_out <- PCA_gen(formated,Pop_ID_Sum)
head(pca_out)
pve <- pca_out$pve[1:5]
pve

pca_df <- pca_out$pca_df[,1:5]
pca_df <- cbind(Pop_ID_Sum,pca_df)

write.csv(pca_df,'~/Documents/Paper1/data_files/pca_df_final_546.csv',row.names = FALSE)
pca_df <- read.csv("~/Documents/Paper1/data_files/pca_df_final_546.csv")



####### PLOT #######

pca_df_B <- merge(ID_transect, pca_df, by="ID")

col_trans <- c("#FFDB6D", "#C4961A","#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

pca_df$Transect <- ID_clim_geo_transect$Transect_SL

pca_B_transect <- ggplot(data = pca_df, aes(x=-1*PC1,y=1*PC2,fill=Transect)) + 
  geom_point(colour='black', size=2, pch=21) + 
  scale_fill_manual(values=colSp2) +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
  theme_bw() + 
  theme(#legend.position = 'none',
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size=13),
    legend.title = element_text(size = 16, colour="black", face="bold", vjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

pca_B_transect
ggsave('~/Documents/poplar_ROAR_files/poplar_files_for_analysis/Paper1_figs/PCA12_transect_labels_546trees_LDpruned_01_10kbwindow_v2.pdf',pca_B_transect,height=4,width=6.5,units='in')





################################## PCA of color-coded parental-types and grey hybrids #########################################

colSp2 <- c("gold", "darkblue", "grey", "deepskyblue")
dadi_pops <- read.table("~/Documents/Paper1/data_files/Pop_ID_3d_dadi_v2.txt")
colnames(dadi_pops) <- c("ID", "Pop")
pca_indv_546_v2 <- right_join(dadi_pops, pca_df, by="ID")
pca_indv_546_v2$Pop = pca_indv_546_v2$Pop.x %>% replace_na('hybrid') 

pca_c <- ggplot(data = pca_indv_546_v2, aes(x=1*PC1,y=-1*PC2,fill=Pop)) + 
  geom_point(colour='black', size=4, pch=21) + 
  scale_fill_manual(values=colSp2) +
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))  +
  theme_bw() + 
  theme(#legend.position = 'none',
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size=13),
    legend.title = element_text(size = 16, colour="black", face="bold", vjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

pca_c
ggsave('~/Documents/Paper1/PCA12_parent_labels_546trees_LDpruned_01_10kbwindow.pdf',pca_c,height=6,width=7,units='in')



pca_4 <- ggplot(data = pca_indv_546_v2, aes(x=1*PC3,y=-1*PC4,fill=Pop)) + 
  geom_point(colour='black', size=4, pch=21) + 
  scale_fill_manual(values=colSp2) +
  xlab(paste("PC",3," (",pve[3]*100,"%)",sep="")) + ylab(paste("PC",4," (",pve[4]*100,"%)",sep=""))  +
  theme_bw() + 
  theme(#legend.position = 'none',
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    legend.text = element_text(size=13),
    legend.title = element_text(size = 16, colour="black", face="bold", vjust=1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

pca_4

################################## PCA of just parental-types #########################################




df012<- fread("~/Documents/Paper1/data_files/poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012",sep="\t",data.table = F)
df012 <- df012[,-1]
df012_indv <- read.table("~/Documents/Paper1/data_files/poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012.indv")
colnames(df012_indv) <- "ID"
dadi_pops <- read.table("~/Documents/Paper1/data_files/Pop_ID_3d_dadi_v2.txt")
colnames(dadi_pops) <- c("ID", "Pop")
dadi_ID <- data.frame(dadi_pops[,1])
colnames(dadi_ID) <- "ID"
df012_labels <- cbind(df012_indv, df012)
df_parents_012 <- inner_join(df012_labels, dadi_ID, by="ID")
df_parents_012 <- df_parents_012[,-1]
head(df_parents_012)

formated<-apply(df_parents_012,2,function(df) as.numeric(df))

Pop_ID_Sum <- dadi_pops
pca_out <- PCA_gen(formated,Pop_ID_Sum)
head(pca_out)
pve <- pca_out$pve[1:5]
pve
#    PC1     PC2     PC3     PC4     PC5 
#.  0.03086 0.01694 0.00885 0.00878 0.00791 

pca_df <- pca_out$pca_df[,1:5]
pca_df <- cbind(Pop_ID_Sum,pca_df)

write.csv(pca_df,'~/Documents/Paper1/pca_df_final_185_dadi_parents.csv',row.names = FALSE)




####### PLOT #######


pca_df <- read.csv('~/Documents/Paper1/pca_df_final_185_dadi_parents.csv')

pve <- c(0.03086, 0.01694, 0.00885, 0.00878, 0.00791 )

pca_parents <- ggplot(data = pca_df, aes(x=-1*PC1,y=PC2, fill=Pop)) + 
  geom_point(size = 8, shape =21) +
  scale_fill_manual(values=c( "gold","darkblue","deepskyblue"))+
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))+
  theme_bw() + 
  theme(#legend.position = 'none',
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    legend.text = element_text(size=13),
    legend.title = element_text(size = 16, colour="black", face="bold", vjust=1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
pca_parents

ggsave('~/Documents/cpDNA/PCA/PCA_dadi_parents_nuclear_LDpruned_snps_v2.pdf',pca_parents,height=6,width=7.5,units='in')








####################################### PCA of cpDNA genetic structure. ########################################

df012<- fread("~/Documents/Paper1/data_files/Populus_plastid_RenameID_185_Daddy_q30_missing05.recode.vcf.012",sep="\t",data.table = F)
df012 <- df012[,-1]
dim(df012). #185 2848

df012_indv <- read.table("~/Documents/Paper1/data_files/Populus_plastid_RenameID_185_Daddy_q30_missing05.recode.vcf.012.indv")
df012_pos <- read.table("~/Documents/Paper1/data_files/Populus_plastid_RenameID_185_Daddy_q30_missing05.recode.vcf.012.pos")
View(df012_indv)


df_gen <- apply(df012, 2, function(df) gsub(-1,NA,df,fixed=TRUE))
View(df_gen)
formated<-apply(df_gen,2,function(df) as.numeric(df))
dim(formated)
View(formated)


Pop_ID_Sum <- dadi_pops
pca_out <- PCA_gen(formated,Pop_ID_Sum)
head(pca_out)
pve <- pca_out$pve[1:5]
pve
#    PC1     PC2     PC3     PC4     PC5 
# 0.23931 0.12291 0.06165 0.03004 0.02335 

pca_df <- pca_out$pca_df[,1:5]
pca_df <- cbind(Pop_ID_Sum,pca_df)

write.csv(pca_df,'~/Documents/Paper1/pca_df_cpDNA_185_dadi_parents.csv',row.names = FALSE)


####### PLOT #######

pve <- c(0.23931, 0.12291, 0.06165, 0.03004, 0.02335 )

pca_cp <- ggplot(data = pca_df, aes(x=-1*PC1,y=PC2, fill=Pop)) + 
  geom_point(size = 8, shape =21) +
  scale_fill_manual(values=c( "gold","darkblue","deepskyblue"))+
  xlab(paste("PC",1," (",pve[1]*100,"%)",sep="")) + ylab(paste("PC",2," (",pve[2]*100,"%)",sep=""))+
  theme_bw() + 
  theme(#legend.position = 'none',
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    legend.text = element_text(size=13),
    legend.title = element_text(size = 16, colour="black", face="bold", vjust=1),
    panel.border = element_rect(size = 1.5, colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())
pca_cp

ggsave('~/Documents/cpDNA/PCA/PCA_dadi_parents_cpDNA_012file.pdf',pca_cp,height=6,width=7.5,units='in')









########### Habitat differences between K=3 parental types #################
install.packages("dunn.test")
library(dunn.test)
install.packages("yarrr")
library(yarrr)
library(dplyr)


data_all <- read.csv("~/Documents/Paper1/data_files/climateNA_7var_1961-1990Y_30arcsec_original_elev.csv")
head(data_all)

dadi_pops <- read.table("~/Documents/Paper1/data_files/Pop_ID_3d_dadi_v2.txt")
head(dadi_pops)
colnames(dadi_pops) <- c("ID", "Pop")

data_dadi_pops <- inner_join(data_all, dadi_pops, by="ID")

data_dadi_balsam <- subset(data_dadi_pops, data_dadi_pops$Pop=="dbalsam")
data_dadi_balsam$Pop_v2 <- "balsam"
data_dadi_coast <- subset(data_dadi_pops, data_dadi_pops$Pop=="coastal_tricho")
data_dadi_coast$Pop_v2 <- "coastal_tricho"
data_dadi_int <- subset(data_dadi_pops, data_dadi_pops$Pop=="interior_tricho")
data_dadi_int$Pop_v2 <- "interior_tricho"
data_dadi_all <- rbind(data_dadi_balsam, data_dadi_coast, data_dadi_int)

###plot elev diff
hist(data_dadi_balsam$elev, breaks=40, col=yarrr::transparent("darkblue",trans.val=0.5), main="Parental-type associations with Elevation", xlab="Elevation (meters)", xlim=c(0,3000))
hist(data_dadi_coast$elev, breaks = 20, add=TRUE, col=yarrr::transparent("gold",0.5))
hist(data_dadi_int$elev, breaks = 20, add=TRUE, col=yarrr::transparent("deepskyblue",0.5))


p <- dunn.test(data_dadi_all$elev, data_dadi_all$Pop, method='bonferroni')


### plot TD diff
hist(data_dadi_balsam$TD, breaks=40, col=yarrr::transparent("darkblue",trans.val=0.5),main="Parental-type associations with TD", xlab="Continentality", xlim=c(10,50))
hist(data_dadi_coast$TD, breaks = 30, add=TRUE, col=yarrr::transparent("gold",0.5))
hist(data_dadi_int$TD, breaks = 10, add=TRUE, col=yarrr::transparent("deepskyblue",0.5))


o <- dunn.test(data_dadi_all$TD, data_dadi_all$Pop, method='bonferroni')


### plot CMD diff
hist(data_dadi_balsam$CMD, breaks=30, col=yarrr::transparent("darkblue",trans.val=0.5), main="Parental-type associations with CMD", xlim=c(0,600), xlab="Climate Moisture Deficit")
hist(data_dadi_coast$CMD, breaks = 20, add=TRUE, col=yarrr::transparent("gold",0.5))
hist(data_dadi_int$CMD, breaks = 30, add=TRUE, col=yarrr::transparent("deepskyblue",0.5))


n <- dunn.test(data_dadi_all$CMD, data_dadi_all$Pop, method='bonferroni')
q <- dunn.test(data_dadi_all$RH, data_dadi_all$Pop, method='bonferroni')
r <- dunn.test(data_dadi_all$lat, data_dadi_all$Pop, method='bonferroni')
s <- dunn.test(data_dadi_all$MAT, data_dadi_all$Pop, method='bonferroni')
t <- dunn.test(data_dadi_all$MAP, data_dadi_all$Pop, method='bonferroni')
u <- dunn.test(data_dadi_all$PAS, data_dadi_all$Pop, method='bonferroni')
v <- dunn.test(data_dadi_all$long, data_dadi_all$Pop, method='bonferroni')

sink('~/Documents/Paper1/significance_test_6climate_3geo_assoc_with_parents_K3_dadi.txt')
o
n
p
q
r
s
t
u
v

sink()


par(mfrow=c(3,3))
boxplot(TD~Pop_v2,data=data_dadi_all, ylab="Continentality", xlab=NULL, ylim=c(0,50))
boxplot(MAT~Pop_v2,data=data_dadi_all, ylab="Mean Annual Temperature", xlab=NULL, ylim=c(-10,15))
boxplot(RH~Pop_v2,data=data_dadi_all, ylab="Relative Humidity", xlab=NULL, ylim=c(40,80))
boxplot(CMD~Pop_v2,data=data_dadi_all,ylab="Climate Moisture Deficit",xlab=NULL, ylim=c(0,800))
boxplot(MAP~Pop_v2,data=data_dadi_all,ylab="Mean Annual Precip.", xlab=NULL, ylim=c(0,2500))
boxplot(PAS~Pop_v2,data=data_dadi_all,ylab="Precip. as snow", xlab=NULL, ylim=c(0,600))
boxplot(elev~Pop_v2,data=data_dadi_all,ylab="Elevation", xlab=NULL, ylim=c(0,3000))
boxplot(lat~Pop_v2,data=data_dadi_all, ylab="Latitude", xlab=NULL, ylim=c(30,70))
boxplot(long~Pop_v2,data=data_dadi_all, ylab="Longitude", xlab=NULL, ylim=c(-160, -90))
dev.copy2pdf(file="~/Documents/Paper1/Boxplots_dadi_parents_K3_distribution_6climate_3geo_vars.pdf", useDingbats=FALSE, family="sans")


par(mfrow=c(3,1))
boxplot(TD~Pop,data=data_dadi_all, ylab="Continentality", ylim=c(0,50))
boxplot(CMD~Pop,data=data_dadi_all,ylab="Climate Moisture Deficit", ylim=c(0,800))
boxplot(elev~Pop,data=data_dadi_all,ylab="Elevation", ylim=c(0,3000))
dev.copy2pdf(file="~/Documents/Paper1/Boxplots_dadi_parents_K3_distribution_climate_assoc.pdf", useDingbats=FALSE, family="sans")

par(mfrow=c(3,1))
boxplot(MAT~Pop,data=data_dadi_all, ylab="Mean Annual Temperature", ylim=c(-10,15))
boxplot(RH~Pop,data=data_dadi_all, ylab="Relative Humidity", ylim=c(40,80))
boxplot(lat~Pop,data=data_dadi_all, ylab="Latitude", ylim=c(30,70))
dev.copy2pdf(file="~/Documents/Paper1/Boxplots_dadi_parents_K3_distribution_climate_assoc_3other_variables.pdf", useDingbats=FALSE, family="sans")
