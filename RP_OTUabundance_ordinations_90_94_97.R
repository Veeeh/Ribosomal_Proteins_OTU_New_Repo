RP.L290C<-read.delim("RP_Cluster_90/Ribosomal_Proteins_S1.1/clustered_otu_out.txt")
RP.L294C<-read.delim("RP_Cluster_94/Ribosomal_Proteins_S1.1/clustered_otu_out.txt")
RP.L297C<-read.delim("RP_Cluster_97/Ribosomal_Proteins_S1.1/clustered_otu_out.txt")
head(RP.L294C)
#RP.L290C_2<-read.delim("RP_Cluster_90/Ribosomal_Proteins_S1.1/clustered_otu_out.txt") testing setseed
colnames(RP.L294C)<-c("otu-num","X06.071", "X06.080", "X06.095", "X06.102", "X06.104", "X06.126", "X06.127",
                   "X06.129", "X06.131", "X06.132", "X06.136", "X06.140", "X06.141", "X06.142", "X06.156",
                   "X06.161", "X06.167", "X06.174", "X06.196", "X06.198", "X06.199", "X06.217", "X06.220",
                   "X07.001", "X07.002", "X07.003", "X07.006", "X07.007", "X07.008", "X07.009", "X07.010",
                   "X07.011", "X07.014", "X07.015", "X07.021", "X07.022", "X07.026", "X07.028", "X07.029",
                   "X07.030", "X07.031", "X07.032", "X07.033", "X07.040", "X07.047", "X07.049", "X07.051",
                   "X07.055", "X07.057", "X08.097", "X08.120", "X08.134", "X08.135", "X08.138", "X08.144",
                   "X08.152", "X08.155", "X08.157", "X08.160", "X08.164", "X08.168", "X08.179", "X08.180",
                   "X08.183", "X08.184", "X08.186","X08.188", "X08.189", "X08.193", "X08.194", "X08.197",
                   "X08.210", "X08.212", "X08.216", "X08.219", "X17.018", "X17.027", "X17.043", "X17.050",
                   "X17.054", "X17.056", "X17.058", "X17.060", "X17.061", "X17.063", "X17.065", "X17.067",
                   "X17.072", "X17.073", "X17.074", "X17.077", "X17.086", "X17.090", "X17.091", "X17.093",
                   "X17.098", "X17.107", "X17.112", "X17.113", "X17.116","seq")
#####remove unwanted column
RP.L290C[1]<-NULL
RP.L294C[1]<-NULL
RP.L297C[1]<-NULL
RP.L290C_2[1]<-NULL##testing set seed!!
######
##Write the files
write.csv(RP.L290C,"C:/Users/V_ONANA/Documents/Ribosomal_Proteins_OTU/RPL2.OTU@90C.csv")
write.csv(RP.L294C,"C:/Users/V_ONANA/Documents/Ribosomal_Proteins_OTU/RPL2.OTU@94C.csv")
write.csv(RP.L297C,"C:/Users/V_ONANA/Documents/Ribosomal_Proteins_OTU/RPL2.OTU@97C.csv")
#write.csv(RP.L290C_2,"C:/Users/V_ONANA/Documents/Ribosomal_Proteins_OTU/RPL2.OTU@90C_NEW_B_each_metagnom.csv") testing set seed 

#Normalise to get relative abundance 
#we normalise using percentage abundance x/sum(x)*100
To_norm.L290<-RP.L290C
To_norm.L290[] <- lapply(To_norm.L290[], function(x){
  # Check if the column is numeric
  if (is.numeric(x)){
    return(x/sum(x))
  } else{
    return(x)
  }
})

############# (NMDS with HI Class and Ecozone)
library(vegan)
library(ggplot2)
library(gplots)
library(vegan)
#make seq rownames for all matrix and change class to matrix 
row.names(To_norm.L290)<-To_norm.L290$seq
To_norm.L290[101]<-NULL
To_norm.L290<-as.matrix(To_norm.L290)
To_norm.L294<-as.matrix(To_norm.L294)
To_norm.L297<-as.matrix(To_norm.L297)
View(To_norm.L290)
is.numeric(To_norm.L290)
#####merge ecozones and HI class
ECO_HI.RP<-merge(Ecozones, HI_class, by = "Lake_ID", all.x = T, all.y = T)
head(ECO_HI.RP)
########## NMDS of RPL290,94,97 ##########
library(vegan)
library(ggplot2)
library(dplyr)
Taxa_RP.NM90<-To_norm.L290
Taxa.NM90.hel<-decostand(t(Taxa_RP.NM90), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM90.NMall<-metaMDS(Taxa.NM90.hel,k=2,trymax=100)
Taxa_NM90.NMall$points


## Building data structure for the graph
Taxa.NM90pos<-data.frame(Taxa_NM90.NMall$points)
View(Taxa.NM90pos)
Taxa.NM90pos2<-cbind(Taxa.NM90pos,rownames(Taxa.NM90pos))
colnames(Taxa.NM90pos2)[3]<-"Lake_ID"
head(Taxa.NM90pos2)
Taxa.NM90graph<-left_join(Taxa.NM90pos2,ECO_HI.RP,by="Lake_ID")
View(Taxa.NM90graph)


########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM90graph, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0, size=2.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-2.0,1.5)+
  ylim(-1.0,1.5)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")

############@ 94 percent 
Taxa_RP.NM94<-To_norm.L294
Taxa.NM94.hel<-decostand(t(Taxa_RP.NM94), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM94.NMall<-metaMDS(Taxa.NM94.hel,k=2,trymax=100)
Taxa_NM94.NMall$points

## Building data structure for the graph
Taxa.NM94pos<-data.frame(Taxa_NM94.NMall$points)
View(Taxa.NM94pos)
Taxa.NM94pos2<-cbind(Taxa.NM94pos,rownames(Taxa.NM94pos))
colnames(Taxa.NM94pos2)[3]<-"Lake_ID"
head(Taxa.NM94pos2)
Taxa.NM94graph<-left_join(Taxa.NM94pos2,ECO_HI.RP,by="Lake_ID")
head(Taxa.NM94graph)

########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM94graph, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0, size=2.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-2.3,1.3)+
  ylim(-1.0,1.2)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")

##############@97 PERECENT 
Taxa_RP.NM97<-To_norm.L297
Taxa.NM97.hel<-decostand(t(Taxa_RP.NM97), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM97.NMall<-metaMDS(Taxa.NM97.hel,k=2,trymax=100)
Taxa_NM97.NMall$points

## Building data structure for the graph
Taxa.NM97pos<-data.frame(Taxa_NM97.NMall$points)
View(Taxa.NM97pos)
Taxa.NM97pos2<-cbind(Taxa.NM97pos,rownames(Taxa.NM97pos))
colnames(Taxa.NM97pos2)[3]<-"Lake_ID"
head(Taxa.NM97pos2)
Taxa.NM97graph<-left_join(Taxa.NM97pos2,ECO_HI.RP,by="Lake_ID")
View(Taxa.NM97graph)

########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM97graph, aes(x = MDS1, y = MDS2)) + 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0,size=2.5)+
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-1.3,1.0)+
  ylim(-1.0,1.0)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")



##########use PCA to show Precent variation
#take transformed data by hellingers 

RP97.rda<-rda(Taxa.NM97.hel)
RP.eigper97<-RP97.rda$CA$eig*100/sum(RP97.rda$CA$eig) # this gets me the list of eigen values in percentage of total variance explained
head(RP.eigper97)
library(tidyverse)

## Building the graph
RP97.PCpos<-RP97.rda$CA$u[,1:2]
RP97.PCpos2<-data.frame(cbind(RP97.PCpos,rownames(RP97.PCpos)))
colnames(RP97.PCpos2)[3]<-"Lake_ID"
RP97.PCgraph<-left_join(RP97.PCpos2,Ecozones,by="Lake_ID") 
library(RColorBrewer)
#KO.PCgraph[,1:2]<-KO.PCpos # for some reason, when I c bind earlier, the numerical values are transformed into factors, so I put them back in

ecozone.colors<-colorRampPalette(brewer.pal(9,"RdYlGn"))(nlevels(RP97.PCgraph$Ecozone))
ecozone.colors<-c("blue","red","orange","purple")

# Creating the plot based on PCA 
ggplot(RP94.PCgraph, 
       aes(y=PC2,x=PC1,color=Ecozone))+
  scale_color_manual(values=ecozone.colors)+
  #geom_text(aes(label=Ecozone),hjust=0, vjust=0)+
  xlim(-0.4,0.4)+
  ylim(-0.4,0.4)+  
  scale_x_discrete(breaks = seq(-0.5, 0.5, by = 0.1))+
  scale_y_discrete(breaks = seq(-0.5, 0.5, by = 0.1))+
  geom_point(size=3)+
  theme_bw() + 
  labs(y=paste("PC2(",as.character(signif(RP.eigper[2],digits=3)),"%)"),
       x=paste("PC1(",as.character(signif(RP.eigper[1],digits=3)),"%)"),
       col="Ecozone")+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")+
  theme(legend.text=element_text(face="bold",size=14), 
        legend.title=element_text(face="bold",size=16),
        plot.title=element_blank(),
        axis.title.x=element_text(colour="black",face="bold",size=16), 
        axis.title.y=element_text(colour="black",face="bold",size=16), 
        axis.text.x=element_text(hjust = 1, colour="black",face="bold",size=14),
        axis.text.y=element_text(colour="black",face="bold",size=14),
        axis.line=element_line(colour="black",size=1,linetype="solid"),
        panel.grid.major = element_blank(), # removes the grid on the graph
        panel.grid.minor= element_blank(),
        panel.border = element_blank(), # this line and the following to remove top and right axis
        panel.background = element_blank(),
        axis.ticks = element_line(colour="black",size=0.5,linetype="solid"))






#################
#######Carry out correlations between percent reads for all groups
###i.e. Bac,Archaea, Euk and Cyano between SSU and RP analysis outputs
#####for comparison 
Correl.matrix<-read.csv("Correlations_percent_SSU_RP_Bac_Euk_Arc_Cya.csv")
Correl.matrix<-Correl.matrix[-c(101:104),]
###############################################
install.packages("ggcorrplot")
library(ggcorrplot)
library(corrplot)
library(RColorBrewer)
install.packages("ggpubr")
library(ggpubr)
library(tidyverse)
install.packages("Hmisc")
library(Hmisc)
library(dplyr)
###########################
#Correlations 
row.names(Correl.matrix)<-Correl.matrix$LAKE_ID
Correl.matrix[1]<-NULL

#make rownames point labels
Correl.matrix$name <- rownames(Correl.matrix)

ggscatter(Correl.matrix, x ="X..Eukarya.SSUrRNA.analysis.", y = "X..Eukarya..RPL2.",
          color = "red", shape = 21, size = 3,fill = "red",
          label = "name", repel=F,
          add = "reg.line", conf.int = TRUE, xlab = "Percent Eukarya SSUrRNA", ylab = "Percent Eukarya RPL2",
          #add.params = list(color = "blue", fill = "lightgray"),
          cor.coef = TRUE, cor.method = "pearson",
          title = "Relationship between percent Eukarya reads in both analysis")
 

library(vegan)
###########Use procrustes to compare ordinations at various clustering levels 
pro_90_94 <- procrustes(X = NMDS90, Y = NMDS94, symmetric = FALSE, scores = "species")
pro

################checking for outliers in metagenomes with boxplots 
###This is done with respect to the number of genes found in them, for example using boxplots:
###Create a matrix for boxplot with number of reads for both analysis
Boxplot_mat<-read.csv("Number_of_reads_SSU_RP_to_create_boxplots.csv")
colnames(Boxplot_mat)<-c("Lake_ID","NOR_SSU","NOR_RP")
Boxplot_mat<-Boxplot_mat[-c(101:104),]
boxplot(Boxplot_mat$NOR_SSU)
boxplot(Boxplot_mat$NOR_SSU)$out
boxplot(Boxplot_mat$NOR_SSU)$stats
#######for RP 
boxplot(Boxplot_mat$NOR_RP)
boxplot(Boxplot_mat$NOR_RP)$out
boxplot(Boxplot_mat$NOR_RP)$stats

###############Note about boxplots 
###From boxplot help:
###Value: List with the following components:
###1.stats: a matrix, each column contains the extreme of the 
###lower whisker, the lower hinge, the median, the upper hinge 
####and the extreme of the upper whisker for one group/plot. 
###If all the inputs have the same class attribute, 
##so will this component.
###2.n:a vector with the number of observations in each group.
###3.conf:a matrix where each column contains the lower and upper 
##extremes of the notch.
###4.out:the values of any data points which lie beyond 
##the extremes of the whiskers.
###5.group: a vector of the same length as out whose elements 
##indicate to which group the outlier belongs.
###6.names:a vector of names for the groups.

#########Do an Envfit on NMDS PLOTs at 90, 94 and 97% clusterings
#####
###Read in Metadata file with correct cations all 100 lakes
Metadata_Rp<-read.csv("Metadata_2017_RP.csv")
#NMDS WITH all parameters for 90C, use points from NMDS ordination
NMDSEnvfit.90C<-envfit(Taxa.NM90pos, env = Metadata_Rp, perm = 999, na.rm = T)#standard envfit
NMDSEnvfit.90C$vectors

##At running all variables, correlated variables were significant (ions) (Landuse), 
##so I will break the metadata into categories, Land use, Water Chemistry and climate
Landuse<-Metadata_Rp[,c(1,11,13:18)]
Water_Chem<-Metadata_Rp[,c(1,20,22,37,39,40:44)]
#merge both into one table
Landuse_Waterchem<-merge(Landuse, Water_Chem, by="Lake_ID", all.x = T, all.y = T)
#NMDS with Landuse_waterchem parameters 90C
#these variables are closely correlating on the plot
##I decided to pick urban, TN and DIC but is this the best approach?
###I should rather do a forward selection so I can pick the most important variables
Landuse_Waterchem_sh<-Landuse_Waterchem[c(1,8,10,13)]
##pick variables that were strong in forward selection
###and make a metadata matrix with these from dbrda at the bottom

NMDSEnvfit.LU90C<-envfit(Taxa.NM90pos, env = Envfit_dbrda_sel, perm = 999, na.rm = T)#standard envfit
NMDSEnvfit.LU90C$vectors
Envfit_dbrda_sel<-Metadata_Rp[,c(1, 8,16,22,40)]
#data for plotting NMDS Landuse90C 
#data for envfit arrows
env.scores.NM.LU90C <- as.data.frame(scores(NMDSEnvfit.LU90C, display = "vectors")) #extracts relevant scores from envifit
env.scores.NM.LU90C <- cbind(env.scores.NM.LU90C, env.variables = rownames(env.scores.NM.LU90C)) #and then gives them their names
#ALL VARIABLES were significant and not correlated only HI and AG @ 0.6 LOW!

##########Plot 
ggplot(data = Taxa.NM90graph, aes(y = MDS2, x = MDS1))+ #sets up the plot. brackets around the entire thing to make it draw automatically
  geom_segment(data = env.scores.NM.LU90C,
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               size=0.5,
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + #arrows for envfit.  doubled the length for similarity to the plot() function. NB check ?envfit regarding arrow length if not familiar with lengths
  geom_text(data = env.scores.NM.LU90C, #labels the environmental variable arrows 
            aes(x = MDS1, y = MDS2, label=env.variables),
            size = 4,
            hjust = 0.5)+
  geom_point(data=Taxa.NM90graph,
             aes(y=MDS2,x=MDS1, shape=Ecozone, colour=HI_Class),
             size=3,
             hjust=0.5)+
  geom_text(data=Taxa.NM90graph,
            aes(label=Lake_ID),size=2.5,colour="black")+
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-2.0,1.5)+
  ylim(-1.0,1.5)+
  #these are the species points, made lighter and a specific shape
  #scale_shape_manual(values = c(1,8,19,5))+ #sets the shape of the plot points instead of using whatever ggplot2 automatically provides
  #coord_cartesian(xlim = c(-0.4,0.4))+  ## NB this changes the visible area of the plot only (this is a good thing, apparently). Can also specify ylim. Here in case you want to set xaxis manually.
  theme_bw()+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")+
  theme(legend.text=element_text(face="bold",size=14), 
        legend.title=element_text(face="bold",size=16),
        plot.title=element_blank(),
        axis.title.x=element_text(colour="black",face="bold",size=16), 
        axis.title.y=element_text(colour="black",face="bold",size=16), 
        axis.text.x=element_text(hjust = 1, colour="black",face="bold",size=14),
        axis.text.y=element_text(colour="black",face="bold",size=14),
        axis.line=element_line(colour="black",size=1,linetype="solid"),
        panel.grid.major = element_blank(), # removes the grid on the graph
        panel.grid.minor= element_blank(),
        panel.border = element_blank(), # this line and the following to remove top and right axis
        panel.background = element_blank(),
        axis.ticks = element_line(colour="black",size=0.5,linetype="solid"))

#######All variables in LU and WC models are significant in the env fitting
###I will use forward selection to pick the most important variables
###DBRDA 
#TEST MODEL ON 90C RPL2
library(ggplot2)
rda.matRPL290<-as.data.frame(t(To_norm.L290))
write.csv(rda.matRPL290,"C:/Users/V_ONANA/Documents/Ribosomal_Proteins_OTU/db_rdaRPL2@90C.csv")
rda.matRPL290<-read.csv("db_rdaRPL2@90C.csv")
row.names(rda.matRPL290)<-rda.matRPL290$Lake_ID
rda.matRPL290[1]<-NULL
#row.names(Landuse_Waterchem)<-Landuse_Waterchem$Lake_ID
row.names(Landuse)<-Landuse$Lake_ID
Landuse[1]<-NULL
row.names(Water_Chem)<-Water_Chem$Lake_ID
Water_Chem[1]<-NULL
Landuse_Waterchem[1]<-NULL
General_Cat<-Metadata_Rp[,c(1,7,8,9,11)]
row.names(General_Cat)<-General_Cat$Lake_ID
General_Cat[1]<-NULL
LU_WCdbrdaRPL290<-capscale(rda.matRPL290~.,data = General_Cat,distance = "bray")
anova(LU_WCdbrdaRPL290)#Check if overall model is significant (Yes-0.001***)
ordiR2step(capscale(rda.matRPL290~1,data=General_Cat),scope=formula(LU_WCdbrdaRPL290),direction="forward",pstep=1000)#f.selection picks best explanatory variables

plot(LU_WCdbrdaRPL290)
summary(LU_WCdbrdaRPL290)
screeplot(LU_WCdbrdaRPL290)

#############plot the dbrda 
###########seperate OTU abundance matrices into euk and bacteria 
####create seperate matrices for different clusters 
To_norm.L290.bac<-To_norm.L290[,-c(14,20,21,26,30,31,38,44,49,54,56,61,67,68,73,79,89,93,94,96,98)]
To_norm.L290.Euk<-To_norm.L290[,c(14,20,21,26,30,31,38,44,49,54,56,61,67,68,73,79,89,93,94,96,98)]
To_norm.L294.bac<-To_norm.L294[,-c(14,20,21,26,30,31,38,44,49,54,56,61,67,68,73,79,89,93,94,96,98)]
To_norm.L294.Euk<-To_norm.L294[,c(14,20,21,26,30,31,38,44,49,54,56,61,67,68,73,79,89,93,94,96,98)]
To_norm.L297.bac<-To_norm.L297[,-c(14,20,21,26,30,31,38,44,49,54,56,61,67,68,73,79,89,93,94,96,98)]
To_norm.L297.Euk<-To_norm.L297[,c(14,20,21,26,30,31,38,44,49,54,56,61,67,68,73,79,89,93,94,96,98)]
###########ordinations @90 forBac only matrix
Taxa_RP.NM90.bac<-To_norm.L290.bac
Taxa.NM90.hel.bac<-decostand(t(Taxa_RP.NM90.bac), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM90.NMall.bac<-metaMDS(Taxa.NM90.hel.bac,k=2,trymax=100)
Taxa_NM90.NMall.bac$points
library(dplyr)

## Building data structure for the graph
Taxa.NM90pos.bac<-data.frame(Taxa_NM90.NMall.bac$points)
View(Taxa.NM90pos.bac)
Taxa.NM90pos2.bac<-cbind(Taxa.NM90pos.bac,rownames(Taxa.NM90pos.bac))
colnames(Taxa.NM90pos2.bac)[3]<-"Lake_ID"
head(Taxa.NM90pos2.bac)
Taxa.NM90graph.bac<-left_join(Taxa.NM90pos2.bac,ECO_HI.RP,by="Lake_ID")
View(Taxa.NM90graph.bac)

########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM90graph.bac, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0, size=2.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-2.0,1.5)+
  ylim(-1.0,1.5)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")

###########ordinations @90 for Euk only matrix
Taxa_RP.NM90.Euk<-To_norm.L290.Euk
Taxa.NM90.hel.Euk<-decostand(t(Taxa_RP.NM90.Euk), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM90.NMall.Euk<-metaMDS(Taxa.NM90.hel.Euk,k=2,trymax=100)
Taxa_NM90.NMall.Euk$points


## Building data structure for the graph
Taxa.NM90pos.Euk<-data.frame(Taxa_NM90.NMall.Euk$points)
View(Taxa.NM90pos.Euk)
Taxa.NM90pos2.Euk<-cbind(Taxa.NM90pos.Euk,rownames(Taxa.NM90pos.Euk))
colnames(Taxa.NM90pos2.Euk)[3]<-"Lake_ID"
head(Taxa.NM90pos2.Euk)
Taxa.NM90graph.Euk<-left_join(Taxa.NM90pos2.Euk,ECO_HI.RP,by="Lake_ID")
View(Taxa.NM90graph.Euk)

########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM90graph.Euk, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0, size=2.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-2.0,1.5)+
  ylim(-1.0,1.5)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")

###########ordinations @94 for Bac only matrix
Taxa_RP.NM94.bac<-To_norm.L294.bac
Taxa.NM94.hel.bac<-decostand(t(Taxa_RP.NM94.bac), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM94.NMall.bac<-metaMDS(Taxa.NM94.hel.bac,k=2,trymax=100)
Taxa_NM94.NMall.bac$points


## Building data structure for the graph
Taxa.NM94pos.bac<-data.frame(Taxa_NM94.NMall.bac$points)
head(Taxa.NM94pos.bac)
Taxa.NM94pos2.bac<-cbind(Taxa.NM94pos.bac,rownames(Taxa.NM94pos.bac))
colnames(Taxa.NM94pos2.bac)[3]<-"Lake_ID"
head(Taxa.NM94pos2.bac)
Taxa.NM94graph.bac<-left_join(Taxa.NM94pos2.bac,ECO_HI.RP,by="Lake_ID")
head(Taxa.NM94graph.bac)

########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM94graph.bac, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0, size=2.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-2.0,1.5)+
  ylim(-1.0,1.5)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")

###########ordinations @94 for Euk only matrix
Taxa_RP.NM94.Euk<-To_norm.L294.Euk
Taxa.NM94.hel.Euk<-decostand(t(Taxa_RP.NM94.Euk), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM94.NMall.Euk<-metaMDS(Taxa.NM94.hel.Euk,k=2,trymax=100)
Taxa_NM94.NMall.Euk$points


## Building data structure for the graph
Taxa.NM94pos.Euk<-data.frame(Taxa_NM94.NMall.Euk$points)
head(Taxa.NM94pos.Euk)
Taxa.NM94pos2.Euk<-cbind(Taxa.NM94pos.Euk,rownames(Taxa.NM94pos.Euk))
colnames(Taxa.NM94pos2.Euk)[3]<-"Lake_ID"
head(Taxa.NM94pos2.Euk)
Taxa.NM94graph.Euk<-left_join(Taxa.NM94pos2.Euk,ECO_HI.RP,by="Lake_ID")
(Taxa.NM94graph.Euk)

########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM94graph.Euk, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0, size=2.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-2.0,2.0)+
  ylim(-0.7,1.6)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")

###########ordinations @97 for Bac only matrix
Taxa_RP.NM97.bac<-To_norm.L297.bac
Taxa.NM97.hel.bac<-decostand(t(Taxa_RP.NM97.bac), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM97.NMall.bac<-metaMDS(Taxa.NM97.hel.bac,k=2,trymax=100)
Taxa_NM97.NMall.bac$points


## Building data structure for the graph
Taxa.NM97pos.bac<-data.frame(Taxa_NM97.NMall.bac$points)
head(Taxa.NM97pos.bac)
Taxa.NM97pos2.bac<-cbind(Taxa.NM97pos.bac,rownames(Taxa.NM97pos.bac))
colnames(Taxa.NM97pos2.bac)[3]<-"Lake_ID"
head(Taxa.NM97pos2.bac)
Taxa.NM97graph.bac<-left_join(Taxa.NM97pos2.bac,ECO_HI.RP,by="Lake_ID")
head(Taxa.NM97graph.bac)

########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM97graph.bac, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0, size=2.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-2.0,1.5)+
  ylim(-1.0,1.5)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")

###########ordinations @97 for Euk only matrix
Taxa_RP.NM97.Euk<-To_norm.L297.Euk
Taxa.NM97.hel.Euk<-decostand(t(Taxa_RP.NM97.Euk), method="hellinger") # I use Hellinger transformation as I have come double zeros, even if not that much!!
Taxa_NM97.NMall.Euk<-metaMDS(Taxa.NM97.hel.Euk,k=2,trymax=100)
Taxa_NM97.NMall.Euk$points


## Building data structure for the graph
Taxa.NM97pos.Euk<-data.frame(Taxa_NM97.NMall.Euk$points)
head(Taxa.NM97pos.Euk)
Taxa.NM97pos2.Euk<-cbind(Taxa.NM97pos.Euk,rownames(Taxa.NM97pos.Euk))
colnames(Taxa.NM97pos2.Euk)[3]<-"Lake_ID"
head(Taxa.NM97pos2.Euk)
Taxa.NM97graph.Euk<-left_join(Taxa.NM97pos2.Euk,ECO_HI.RP,by="Lake_ID")
View(Taxa.NM97graph.Euk)

########Now make both HI Class and Eco in one graph
ggplot(Taxa.NM97graph.Euk, aes(x = MDS1, y = MDS2)) + 
  geom_point(size = 3, aes( shape =Ecozone, colour = HI_Class))+ 
  geom_text(aes(label=Lake_ID),hjust=0, vjust=0, size=2.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "HI Class", y = "NMDS2", shape = "Ecozone")  + 
  scale_colour_manual(values = c("red", "green","blue"))+
  xlim(-1.7,1.7)+
  ylim(-0.7,1.5)+
  geom_vline(xintercept = 0,linetype="dotted")+
  geom_hline(yintercept = 0,linetype="dotted")
memory.limit()
###increase storage capacity from 16227 to 3 times of it 
memory.limit(size=48681)
version
###########return capacity 
memory.limit(size = 16227) #could not decrese memory limit!!!!!!!

############Now, I have identified that random seed changes OTU numbers so I set the seed
#########to 99999 to have a constant number accross all matrices (Taxa, OTU)
###########to allow comparison 
########read in the OTU at 90C again 
RPL2.90C.set.seed<-read.delim("RP_Cluster_90/Ribosomal_Proteins_S1.1/clustered_otu_out.txt")
#####now number is same for OTU and Taxa table 
######Change col names 
colnames(RPL2.90C.set.seed)<-c("otu-num","X06.071", "X06.080", "X06.095", "X06.102", "X06.104", "X06.126", "X06.127",
                      "X06.129", "X06.131", "X06.132", "X06.136", "X06.140", "X06.141", "X06.142", "X06.156",
                      "X06.161", "X06.167", "X06.174", "X06.196", "X06.198", "X06.199", "X06.217", "X06.220",
                      "X07.001", "X07.002", "X07.003", "X07.006", "X07.007", "X07.008", "X07.009", "X07.010",
                      "X07.011", "X07.014", "X07.015", "X07.021", "X07.022", "X07.026", "X07.028", "X07.029",
                      "X07.030", "X07.031", "X07.032", "X07.033", "X07.040", "X07.047", "X07.049", "X07.051",
                      "X07.055", "X07.057", "X08.097", "X08.120", "X08.134", "X08.135", "X08.138", "X08.144",
                      "X08.152", "X08.155", "X08.157", "X08.160", "X08.164", "X08.168", "X08.179", "X08.180",
                      "X08.183", "X08.184", "X08.186","X08.188", "X08.189", "X08.193", "X08.194", "X08.197",
                      "X08.210", "X08.212", "X08.216", "X08.219", "X17.018", "X17.027", "X17.043", "X17.050",
                      "X17.054", "X17.056", "X17.058", "X17.060", "X17.061", "X17.063", "X17.065", "X17.067",
                      "X17.072", "X17.073", "X17.074", "X17.077", "X17.086", "X17.090", "X17.091", "X17.093",
                      "X17.098", "X17.107", "X17.112", "X17.113", "X17.116","seq")
#####write file 
write.csv(RPL2.90C.set.seed,"C:/Users/V_ONANA/Documents/Ribosomal_Proteins_OTU/RP_Cluster_90/RPL2.OTU@90C.set.seed.csv")
