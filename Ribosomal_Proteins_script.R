RP.S1<-read.delim("Ribosomal_Proteins_S1.1/otu_out.txt")

RP.S1.c<-read.delim("Ribosomal_Proteins_S1.1/clustered_otu_out.txt")
head(RP.S1)
colnames(RP.S1)<-c("seq", "X06.071", "X06.080", "X06.095", "X06.102", "X06.104", "X06.126", "X06.127",
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
 "X17.098", "X17.107", "X17.112", "X17.113", "X17.116")
names.test<-as.data.frame(names(RP.S1))
library(tidyverse)
library(dplyr)
library(tidyr)
test.2<-separate(p,names.RP.S1., into=c("L","R"), sep="_")
p <- transform(names.test, class=as.numeric(as.character(names.test)))
test.2[c(2,3)]<-NULL
test.2$L

#Normalise to get relative abundance 
#subset the data but before then, we normalise using percentage abundance x/sum(x)*100
To_norm.RP<-RP.S1
To_norm.RP[] <- lapply(To_norm.RP[], function(x){
  # Check if the column is numeric
  if (is.numeric(x)){
    return(x/sum(x))
  } else{
    return(x)
  }
})

#############oRDInations
library(vegan)
library(ggplot2)
library(gplots)
row.names(To_norm.RP)<-To_norm.RP$seq
To_norm.RP[1]<-NULL

#transform data by hellingers 
RP.hel<-decostand(t(To_norm.RP),method="hellinger") 
RP.rda<-rda(RP.hel)
RP.eigper<-RP.rda$CA$eig*100/sum(RP.rda$CA$eig) # this gets me the list of eigen values in percentage of total variance explained
biplot(RP.rda)
Ecozones<-read.csv("Ecozones_Ribosomal_P.csv") 
## Building the graph
RP.PCpos<-RP.rda$CA$u[,1:2]
RP.PCpos2<-data.frame(cbind(RP.PCpos,rownames(RP.PCpos)))
colnames(RP.PCpos2)[3]<-"Lake_ID"
RP.PCgraph<-left_join(RP.PCpos2,Ecozones,by="Lake_ID") 
library(RColorBrewer)
#KO.PCgraph[,1:2]<-KO.PCpos # for some reason, when I c bind earlier, the numerical values are transformed into factors, so I put them back in

ecozone.colors<-colorRampPalette(brewer.pal(9,"RdYlGn"))(nlevels(RP.PCgraph$Ecozone))
ecozone.colors<-c("blue","red","orange","purple")

# Creating the plot based on PCA 
ggplot(RP.PCgraph, 
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

####use Human impact class for ordination
## Building the graph
HI_class<-read.csv("HI_Class_recategorised_RP.csv")
RP.PCpos<-RP.rda$CA$u[,1:2]
RP.PCpos2<-data.frame(cbind(RP.PCpos,rownames(RP.PCpos)))
colnames(RP.PCpos2)[3]<-"Lake_ID"
RP.PCpos.HI<-RP.PCpos2
RP.PCgraph.HI<-left_join(RP.PCpos.HI,HI_class,by="Lake_ID") 
library(RColorBrewer)
#KO.PCgraph[,1:2]<-KO.PCpos # for some reason, when I c bind earlier, the numerical values are transformed into factors, so I put them back in

ecozone.colors.HI<-colorRampPalette(brewer.pal(9,"RdYlGn"))(nlevels(RP.PCgraph.HI$HI_Class))
ecozone.colors.HI<-c("red","green","blue")

# Creating the plot based on PCA 
ggplot(RP.PCgraph.HI, 
       aes(y=PC2,x=PC1,color=HI_Class))+
  scale_color_manual(values=ecozone.colors.HI)+
  #geom_text(aes(label=Ecozone),hjust=0, vjust=0)+
  xlim(-0.4,0.4)+
  ylim(-0.4,0.4)+  
  scale_x_discrete(breaks = seq(-0.5, 0.5, by = 0.1))+
  scale_y_discrete(breaks = seq(-0.5, 0.5, by = 0.1))+
  geom_point(size=3)+
  theme_bw() + 
  labs(y=paste("PC2(",as.character(signif(RP.eigper[2],digits=3)),"%)"),
       x=paste("PC1(",as.character(signif(RP.eigper[1],digits=3)),"%)"),
       col="HI_Class")+
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

#RP.S2
RP.S2<-read.delim("Ribosomal_Proteins_S1.2/otu_out.txt")

RP.S3<-read.delim("Ribosomal_Proteins_S1.3/otu_out.txt")
head(RP.S2)
colnames(RP.S2)<-c("seq", "X06.071", "X06.080", "X06.095", "X06.102", "X06.104", "X06.126", "X06.127",
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
                   "X17.098", "X17.107", "X17.112", "X17.113", "X17.116")


#Normalise to get relative abundance 
#subset the data but before then, we normalise using percentage abundance x/sum(x)*100
To_norm.RP.2<-RP.S2
To_norm.RP.2[] <- lapply(To_norm.RP.2[], function(x){
  # Check if the column is numeric
  if (is.numeric(x)){
    return(x/sum(x))
  } else{
    return(x)
  }
})
row.names(To_norm.RP.2)<-To_norm.RP.2$seq
To_norm.RP.2[1]<-NULL
#transform data by hellingers 
RP.hel.S2<-decostand(t(To_norm.RP.2),method="hellinger") 
RP.rda.2<-rda(RP.hel.S2)
RP.eigper.2<-RP.rda.2$CA$eig*100/sum(RP.rda.2$CA$eig) # this gets me the list of eigen values in percentage of total variance explained
biplot(RP.rda.2)
Ecozones<-read.csv("Ecozones_Ribosomal_P.csv") 
## Building the graph
RP.S2.PCpos<-RP.rda.2$CA$u[,1:2]
RP.S2.PCpos2<-data.frame(cbind(RP.S2.PCpos,rownames(RP.S2.PCpos)))
colnames(RP.S2.PCpos2)[3]<-"Lake_ID"
RP.S2.PCgraph<-left_join(RP.S2.PCpos2,Ecozones,by="Lake_ID") 
library(RColorBrewer)
#KO.PCgraph[,1:2]<-KO.PCpos # for some reason, when I c bind earlier, the numerical values are transformed into factors, so I put them back in

ecozone.colors<-colorRampPalette(brewer.pal(9,"RdYlGn"))(nlevels(RP.PCgraph$Ecozone))
ecozone.colors<-c("blue","red","orange","purple")

# Creating the plot based on PCA 
ggplot(RP.S2.PCgraph, 
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

####use Human impact class for ordination
## Building the graph
HI_class<-read.csv("HI_Class_recategorised_RP.csv")
RP.PCpos<-RP.rda$CA$u[,1:2]
RP.PCpos2<-data.frame(cbind(RP.PCpos,rownames(RP.PCpos)))
colnames(RP.PCpos2)[3]<-"Lake_ID"
RP.PCpos.HI<-RP.PCpos2
RP.PCgraph.HI<-left_join(RP.PCpos.HI,HI_class,by="Lake_ID") 
library(RColorBrewer)
#KO.PCgraph[,1:2]<-KO.PCpos # for some reason, when I c bind earlier, the numerical values are transformed into factors, so I put them back in

ecozone.colors.HI<-colorRampPalette(brewer.pal(9,"RdYlGn"))(nlevels(RP.PCgraph.HI$HI_Class))
ecozone.colors.HI<-c("red","green","blue")

# Creating the plot based on PCA 
ggplot(RP.PCgraph.HI, 
       aes(y=PC2,x=PC1,color=HI_Class))+
  scale_color_manual(values=ecozone.colors.HI)+
  #geom_text(aes(label=Ecozone),hjust=0, vjust=0)+
  xlim(-0.4,0.4)+
  ylim(-0.4,0.4)+  
  scale_x_discrete(breaks = seq(-0.5, 0.5, by = 0.1))+
  scale_y_discrete(breaks = seq(-0.5, 0.5, by = 0.1))+
  geom_point(size=3)+
  theme_bw() + 
  labs(y=paste("PC2(",as.character(signif(RP.eigper[2],digits=3)),"%)"),
       x=paste("PC1(",as.character(signif(RP.eigper[1],digits=3)),"%)"),
       col="HI_Class")+
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

##################################################
#################################

RP.S3<-read.delim("Ribosomal_Proteins_S1.3/otu_out.txt")
RP.S4<-read.delim("Ribosomal_Proteins_S1.4/otu_out.txt")
RP.S5<-read.delim("Ribosomal_Proteins_S1.5/otu_out.txt")
RP.S6<-read.delim("Ribosomal_Proteins_S1.6/otu_out.txt")
RP.S7<-read.delim("Ribosomal_Proteins_S1.7/otu_out.txt")
RP.S8<-read.delim("Ribosomal_Proteins_S1.8/otu_out.txt")
RP.S9<-read.delim("Ribosomal_Proteins_S1.9/otu_out.txt")
RP.S10<-read.delim("Ribosomal_Proteins_S1.10/otu_out.txt")
RP.S11<-read.delim("Ribosomal_Proteins_S1.11/otu_out.txt")
RP.S12<-read.delim("Ribosomal_Proteins_S1.12/otu_out.txt")
RP.S13<-read.delim("Ribosomal_Proteins_S1.13/otu_out.txt")
RP.S14<-read.delim("Ribosomal_Proteins_S1.14/otu_out.txt")

RP.S3.c<-read.delim("Ribosomal_Proteins_S1.3/clustered_otu_out.txt")
dim(RP.S1)
dim(RP.S2)
dim(RP.S4)
dim(RP.S5)
dim(RP.S6)
dim(RP.S7)
dim(RP.S8)
dim(RP.S9)
dim(RP.S10)
dim(RP.S11)
dim(RP.S12)
dim(RP.S13)
dim(RP.S14)
colnames(RP.S13)<-c("seq", "X06.071", "X06.080", "X06.095", "X06.102", "X06.104", "X06.126", "X06.127",
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
                   "X17.098", "X17.107", "X17.112", "X17.113", "X17.116")



