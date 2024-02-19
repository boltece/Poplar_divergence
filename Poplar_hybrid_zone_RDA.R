
#### Bolte et al., 2024
#### Populus trichocarpa x balsamifera hybrid zone 

## Dr. Connie Bolte; conniebolte@gmail.com




######################################## RDA ########################################

#straight from forester paper
#https://popgen.nescent.org/2018-03-27_RDA_GEA.html
install.packages("devtools")
library(devtools)
install.packages("data.table")
library(data.table)
install.packages("tidyverse")
library(tidyverse)
install.packages("vegan")
library(vegan)
install.packages("psych")
library(psych)
library(viridis)
install.packages("ggsci")
library(ggsci)
install.packages("ggcorrplot")
library(ggcorrplot)
require(ggrepel)
library(lme4)
library(readr)


############################### ENV data for RDA ##################################
setwd("./RDA/")

df_mod <- read.csv("data_files/climateNA_7var_1961-1990Y_30arcsec_original_elev.csv")


head(df_mod)

MAT <- df_mod$MAT
TD <- df_mod$TD
MAP <- df_mod$MAP
RH <- df_mod$RH
PAS <- df_mod$PAS
CMD <- df_mod$CMD


climate_6var <- cbind(TD, MAP, MAT, PAS, CMD, RH)
clim6_scaled <- apply(climate_6var, 2, scale)
clim6_df <-as.data.frame(clim6_scaled)

lat <- df_mod$lat
long <- df_mod$long
elev <- df_mod$elev
geo_df <- cbind(long, lat, elev)
geo_scaled <- apply(geo_df,2, scale)
geo_scaled <- as.data.frame(geo_scaled)

geo_clim <- cbind(geo_scaled, clim6_df)


################################ GENETIC data for RDA ##################################


df012<- fread("./poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012",sep="\t",data.table = F)
df012 <- df012[,-1]
dim(df012)

formated<-apply(df012,2,function(df) as.numeric(df))
dim(formated)


### writing genetic data to a text file 
df012_txt <- df012
df012_indv <- read.table("./poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012.indv")

rownames(df012_txt)<- df012_indv$V1
View(df012_txt[1:100,1:100])

df012_pos <- read.table("./poplar_546trees_mac2_missingNONE_LDpruned_POStxtFILE.recode.vcf.012.pos")



############################### SCALE genetic data ####################################

colmean<-apply(formated,2,mean,na.rm=TRUE) 
## all missing data was filtered out using vcftools prior to analysis so na.rm=TRUE is irrelevant for this df

normalize<-matrix(nrow = nrow(formated),ncol=ncol(formated))

af<-colmean/2

for (m in 1:length(af)){
  nr<-formated[ ,m]-colmean[m]
  dn<-sqrt(af[m]*(1-af[m]))
  normalize[ ,m]<-nr/dn
}


scaled_df012 <- as.data.frame(normalize)
View(scaled_df012[1:20,1:20])




########################### partitioning of variance ################################



vp<- varpart(scaled_df012, ~ as.matrix(clim6_df), ~ as.matrix(geo_scaled))
vp
sink('./...vp_analysis-output_6climate_variables_3geo_30sec.txt')
vp
sink()



climate_ind_accounts <- (0.00593/0.02174)*100 # % PVE
climate_ind_accounts #27.27691
geo_ind_accounts <- (0.00305/0.02174)*100 # % PVE
geo_ind_accounts # 14.02944
confounded_effect <- ((0.02174-0.01869-0.01581)/0.02174)*100
confounded_effect # -58.69365




#### RDA of geo and clim ####

n <- rda(formula = scaled_df012 ~.,scale=FALSE, center=TRUE, data = geo_clim)


anova(n)


saveRDS(n,'./RDA_poplar_546trees_6clim_3geo_scaled.RDS')

RsquareAdj(n)

summary(eigenvals(n,model='constrained'))

n <- readRDS("data_files/RDA_poplar_546trees_6clim_3geo_scaled.RDS")

RDA1_gen <- sort(abs(summary(n)$biplot[,1]),decreasing=TRUE)
RDA1_gen
#       TD       MAT       lat       CMD       MAP        RH      elev      long       PAS 
#  0.8733938 0.8361764 0.6858393 0.5986085 0.3219190 0.3044263 0.2222795 0.2039641 0.1232221 

RDA2_gen <- sort(abs(summary(n)$biplot[,2]),decreasing=TRUE)
RDA2_gen
#   RH       elev      long       lat       MAP       CMD       PAS        TD       MAT 
# 0.7626349 0.7340418 0.6985703 0.6381888 0.5080226 0.5064831 0.4519596 0.2507629 0.1743819 

sort(abs(summary(n)$biplot[,1]) + abs(summary(n)$biplot[,2]),decreasing=TRUE)
#       lat        TD       CMD        RH       MAT      elev      long       MAP       PAS 
#.  1.3240281 1.1241567 1.1050916 1.0670611 1.0105583 0.9563213 0.9025344 0.8299417 0.5751817 

RDA_df <- as.data.frame(summary(n)$biplot)
write.csv(RDA_df,'./RDA_poplar_546trees_6clim_3geo_scaled.csv',row.names = F)


anova(n,by="axis")


###################################### ggplot ############################################


n <- readRDS('./RDA_poplar_546trees_6clim_3geo_scaled.RDS')


### file prep 

n_sum_pve <- summary(n)

RDA1_pve <- paste("RDA1 (",round((n_sum_pve$concont$importance[2,1]*100),digits=2),"%)", sep="")
RDA2_pve <- paste("RDA2 (",round((n_sum_pve$concont$importance[2,2]*100),digits=2),"%)", sep="")
RDA3_pve <- paste("RDA3 (",round((n_sum_pve$concont$importance[2,3]*100),digits=2),"%)", sep="")
n_sum <- summary(n)
n_sum <- scores(n, display=c("sp","sites", "bp")) 
rda_snp <- as.data.frame(n_sum$species)
rda_indv <- as.data.frame(n_sum$sites)


pca_df <- read.csv("data_files/pca_df_final_546.csv")

rda_df2 <- data.frame(ID=as.character(df_mod$ID),ANC=as.character(pca_df$Pop),
                      TRANSECT=as.character(df_mod$Transect_SL)) 


              
rda_indv <- cbind(rda_df2,rda_indv)
write.csv(rda_indv, "./RDA_indv_assignment_loadings_6clim_3geo_scaled_modified_transects.csv")
rda_indv <-read.csv("./RDA_indv_assignment_loadings_6clim_3geo_scaled_modified_transects.csv")
rda_biplot <- as.data.frame(n_sum$biplot)
rda_biplot$var <- row.names(rda_biplot)

basplot <- plot(n)
mult <- attributes(basplot$biplot)$arrow.mul





####  LABEL PLOT with parental-types and grey hybrids ####

rda_indv_546 <- read.csv("~/Documents/Paper1/data_files/RDA_indv_assignment_loadings_6clim_3geo_scaled_modified_transects.csv")
head(rda_indv_546)
rda_indv_546 <- rda_indv_546[,-1]
dadi_pops <- read.table("~/Documents/Paper1/data_files/Pop_ID_3d_dadi_v2.txt")
colnames(dadi_pops) <- c("ID", "Pop")
rda_indv_546_v2 <- right_join(dadi_pops, rda_indv_546, by="ID")
rda_indv_546_v2$Pop = rda_indv_546_v2$Pop %>% replace_na('hybrid') 
write.table(rda_indv_546_v2, "~/Documents/Paper1/Pop_IDs_dadi_samples_and_hybrids_labeled.txt", sep="\t", quote=F, col.names = T, row.names = F)
rda_biplot_n<- as.data.frame(n_sum$biplot)
rda_biplot_n$var <- row.names(rda_biplot_n)

basplot <- plot(n)
mult <- attributes(basplot$biplot)$arrow.mul

colSp2 <- c("gold", "darkblue", "grey", "deepskyblue")

rda_plot_n <- ggplot(data=rda_indv_546_v2, aes(x=RDA1, y=RDA2)) + 
  geom_point(data=rda_indv_546_v2, aes(x=RDA1, y=RDA2,fill=Pop),pch=21,col='black',size=3) +
  xlim(-16, 16) +
  ylim(-16,15) +
  scale_fill_manual(values=colSp2) +
  geom_segment(data = rda_biplot_n,
               aes(x = 0, xend = mult * RDA1 * 1.5,y = 0, yend = mult * RDA2 * 1.5),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot_n,
                   aes(x= (mult + mult/2) * RDA1, y = (mult + mult/2) * RDA2, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 4,fontface = "bold") + 
  xlab(RDA1_pve) + ylab(RDA2_pve) +theme_bw()  + 
  theme(#legend.position = "none",
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

rda_plot_n

ggsave('~/Documents/Paper1/RDA_546samples_3D_parents_color_coded_hybrids_grey_clim_geo_scaled_30sec_NOconditioning_v2.pdf',rda_plot_n,height=5,width=7,units='in')
ggsave('~/Documents/Paper1/RDA_546samples_3D_parents_color_coded_hybrids_grey_clim_geo_scaled_30sec_NOconditioning_v3.pdf',rda_plot_n,height=6,width=7,units='in')




### label by transect

col_trans <- c("#FFDB6D", "#C4961A","#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")

rda_plot_transect <- ggplot(data=rda_indv, aes(x=RDA1, y=RDA2)) + 
  geom_point(data=rda_indv, aes(x=RDA1, y= RDA2,fill=TRANSECT),pch=21,col='black',size=3) +
  xlim(-15, 15) +
  ylim(-15,15) +
  scale_fill_manual(values=col_trans) +
  geom_segment(data = rda_biplot,
               aes(x = 0, xend = mult * RDA1 * 0.95,y = 0, yend = mult * RDA2 * 0.95), linewidth=1.0,
               arrow = arrow(length = unit(0.35, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot,
                   aes(x= (mult + mult/7) * RDA1, y = (mult + mult/7) * RDA2, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 5,fontface = "bold") + 
  xlab(RDA1_pve) + ylab(RDA2_pve) +theme_bw()  + 
  theme(#legend.position = "none",
    axis.text = element_text(size=12), 
    axis.title = element_text(size = 14, colour="black",face = "bold",vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

rda_plot_transect 

ggsave('./RDA_Transects_label_V2_6clim_3geo_scaled.pdf',rda_plot_transect,height=6,width=9,units='in')




### label by K = 2 ancestral coefficient


rda_plot <- ggplot(data=rda_indv, aes(x=RDA1, y=RDA2)) + 
  geom_point(data=rda_indv, aes(x=RDA1, y=RDA2,fill=ANC),pch=21,col='black',size=3) +
  scale_fill_gradient(low="deepskyblue", high="blue4") +
  xlim(-15, 15) +
  ylim(-15,15) +
  geom_segment(data = rda_biplot,
               aes(x = 0, xend = mult * RDA1 * 0.95,y = 0, yend = mult * RDA2 * 0.95),
               arrow = arrow(length = unit(0.35, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot,
                   aes(x= (mult + mult/5) * RDA1, y = (mult + mult/5) * RDA2, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 5,fontface = "bold") + 
  xlab(RDA1_pve) + ylab(RDA2_pve) +theme_bw()  + 
  theme(#legend.position = "none",
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

rda_plot  

ggsave('./RDA_label_6clim_3geo_scaled.pdf',rda_plot,height=6,width=9,units='in')









########################## RDA cpDNA ###########################

df012<- fread("~/Documents/Paper1/data_files/Populus_plastid_RenameID_185_Daddy_q30_missing05.recode.vcf.012",sep="\t",data.table = F)
df012 <- df012[,-1]
dim(df012)

df012_indv <- read.table("~/Documents/Paper1/data_files/Populus_plastid_RenameID_185_Daddy_q30_missing05.recode.vcf.012.indv")
df012_pos <- read.table("~/Documents/Paper1/data_files/Populus_plastid_RenameID_185_Daddy_q30_missing05.recode.vcf.012.pos")

df_gen <- apply(df012, 2, function(df) gsub(-1,NA,df,fixed=TRUE))
View(df_gen)
formated<-apply(df_gen,2,function(df) as.numeric(df))
dim(formated)
View(formated)

colmean<-apply(formated,2,mean,na.rm=TRUE)

normalize<-matrix(nrow = nrow(formated),ncol=ncol(formated))

af<-colmean/2

for (m in 1:length(af)){
  nr<-formated[ ,m]-colmean[m]
  dn<-sqrt(af[m]*(1-af[m]))
  normalize[ ,m]<-nr/dn
}

normalize[is.na(normalize)]<-0

scaled_df012 <- as.data.frame(normalize)
View(scaled_df012[1:20,1:20])


##### climate geo  #######
climate_geo_csv <- read.csv("~/Documents/Paper1/data_files/climateNA_7var_1961-1990Y_30arcsec_original_elev.csv")
dadi_pops <- read.table("~/Documents/Paper1/data_files/Pop_ID_3d_dadi_v2.txt")
colnames(dadi_pops) <- c("ID", "Pop")
dadi_ID <- data.frame(dadi_pops[,1])
colnames(dadi_ID) <- "ID"
df_parents_clim <- inner_join(climate_geo_csv, dadi_ID, by="ID")

MAT <- df_parents_clim$MAT
TD <- df_parents_clim$TD
MAP <- df_parents_clim$MAP
RH <- df_parents_clim$RH
PAS <- df_parents_clim$PAS
CMD <- df_parents_clim$CMD
dadi_climate_6var <- cbind(MAT,TD, MAP, PAS, CMD, RH)
dadi_clim6_scaled <- apply(dadi_climate_6var, 2, scale)
clim_df <-as.data.frame(dadi_clim6_scaled)


lat <- df_parents_clim$lat
long <- df_parents_clim$long
elev <- df_parents_clim$elev
dadi_geo_3var <- cbind(lat, long, elev)
dadi_geo3_scaled <- apply(dadi_geo_3var, 2, scale)
geo_df <- as.data.frame(dadi_geo3_scaled)


clim_geo_df <- cbind(clim_df,geo_df)


### just climate RDA
cp <- rda(formula = scaled_df012 ~ MAT + CMD + RH + TD + PAS + MAP + lat + long + elev ,scale=FALSE, center=TRUE, data = clim_geo_df)
cp
sink('~/Documents/Paper1/cpDNA/cpDNA_rda_dadi_3D_parents_NOconditioning_6clim_3geo_variables_30sec.txt')
cp
sink()

o5_sum_pve <- summary(cp)

RDA1_pve5 <- paste("RDA1 (",round((o5_sum_pve$concont$importance[2,1]*100),digits=2),"%)", sep="")
RDA2_pve5 <- paste("RDA2 (",round((o5_sum_pve$concont$importance[2,2]*100),digits=2),"%)", sep="")
o5_sum <- summary(cp)
o5_sum <- scores(cp, display=c("sp","sites", "bp")) 
rda_snp5 <- as.data.frame(o5_sum$species)
rda_indv5 <- as.data.frame(o5_sum$sites)
rda_indv5 <- cbind(dadi_pops,rda_indv5)
write.csv(rda_indv5, "~/Documents/Paper1/cpDNA/cpDNA_RDA_indv_assignment_loadings_dadi_3D_6clim_3geo_NOconditioning_30sec.csv")
rda_biplot5<- as.data.frame(o5_sum$biplot)
rda_biplot5$var <- row.names(rda_biplot5)

basplot5 <- plot(cp)
mult <- attributes(basplot5$biplot)$arrow.mul

colSp <- c("gold", "darkblue", "deepskyblue")

rda_plot_o5 <- ggplot(data=rda_indv5, aes(x=RDA1, y=RDA2)) + 
  geom_point(data=rda_indv5, aes(x=RDA1, y=RDA2,fill=Pop),pch=21,col='black',size=4) +
  xlim(-5, 5) +
  ylim(-9,3) +
  scale_fill_manual(values=colSp) +
  geom_segment(data = rda_biplot5,
               aes(x = 0, xend = mult * RDA1 * 1.5,y = 0, yend = mult * RDA2 * 1.5),
               arrow = arrow(length = unit(0.25, "cm")), colour = "black") + #grid is required for arrow to work.
  geom_label_repel(data = rda_biplot5,
                   aes(x= (mult + mult/2) * RDA1, y = (mult + mult/2) * RDA2, #we add 10% to the text to push it slightly out from arrows
                       label = var), #otherwise you could use hjust and vjust. I prefer this option
                   size = 4,fontface = "bold") + 
  xlab(RDA1_pve5) + ylab(RDA2_pve5) +theme_bw()  + 
  theme(#legend.position = "none",
    axis.text = element_text(size=13), 
    axis.title = element_text(size = 16, colour="black",face = "bold",vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank())

rda_plot_o5

ggsave('~/Documents/Paper1/cpDNA/cpDNA_RDA_dadi_3D_parents_label_clim_geo_scaled_30sec_NOconditioning.pdf',rda_plot_o5,height=6,width=7,units='in')




saveRDS(cp,'~/Documents/Paper1/data_files/cpDNA_RDA_poplar_3D_dadi_parents_6clim_3geo_scaled.RDS')

RsquareAdj(cp)


summary(eigenvals(cp,model='constrained'))



RDA1_gen <- sort(abs(summary(cp)$biplot[,1]),decreasing=TRUE)
RDA1_gen
#      TD        RH       MAP       MAT       PAS       lat      elev      long       CMD 
#.  0.8182202 0.7760410 0.7336191 0.6839999 0.5975862 0.2526448 0.2062524 0.1763108 0.1393823 

RDA2_gen <- sort(abs(summary(cp)$biplot[,2]),decreasing=TRUE)
RDA2_gen
#        CMD       elev        lat         RH        MAT       long         TD        MAP        PAS 
#.  0.92135047 0.71568877 0.64723782 0.56755459 0.33648740 0.31960567 0.22441709 0.18179827 0.05038441 

sort(abs(summary(cp)$biplot[,1]) + abs(summary(cp)$biplot[,2]),decreasing=TRUE)
#.    RH       CMD        TD       MAT      elev       MAP       lat       PAS      long 
#1.3435956 1.0607327 1.0426373 1.0204873 0.9219412 0.9154174 0.8998827 0.6479706 0.4959164 


anova(cp)



#### variance partitioning ####

vp_cp<- varpart(scaled_df012, ~ as.matrix(clim_df), ~ as.matrix(geo_df))
vp_cp

sink('~/Documents/Paper1/cpDNA/cpDNA_VP_dadi_3D_parents_6clim_3geo_output_30sec.txt')
vp_cp
sink()



climate_ind_accounts <- (0.12240/0.30088)*100 # % PVE
climate_ind_accounts # 40.6807
geo_ind_accounts <- (0.00864/0.30088)*100 # % PVE
geo_ind_accounts # 2.8716
confounded_effect <- ((0.30088-0.29224-0.17848)/0.30088)*100
#-56.447753


