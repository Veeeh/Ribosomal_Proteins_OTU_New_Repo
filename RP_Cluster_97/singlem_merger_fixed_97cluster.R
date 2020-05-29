#!/usr/bin/env Rscript


library("optparse")

option_list = list(
  make_option(c("-p", "--path"), type="character", default=NULL, 
              help="dataset path name (should only contain output files for a single protein", metavar="character"),
  make_option(c("-o", "--otu"), type="character", default="otu_out.txt", 
              help="otu file name [default= %default]", metavar="character"),
  make_option(c("-t", "--threshold"), type="character", default="0.97", 
              help="threshold id for clustering [default= %default]", metavar="character"),
  make_option(c("-c", "--cluster"), type="character", default="clustered_otu_out.txt", 
              help="clustered otu file name [default= %default]", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$path)){
  print_help(opt_parser)
  stop("Path name is missing.n", call.=FALSE)
}

###Program
library(ape,quietly=TRUE)
library(kmer,quietly=TRUE)
library(phylotools,quietly=TRUE)

path<-(opt$path)
names<-list.files(as.character(opt$path),pattern=".unifrac")

#Initialize OTU table
df<-read.table(paste0(path,names[1]))
df<-df[c(1)]
colnames(df)<-c("seq")

#Add OTUs from all datasets in the folder
for (i in names) {
  var<-read.table(paste0(path,i))
  var2<-var[c(1,3)]
  colnames(var2)<-c("seq",paste(i))
  #print(head(var2))
  df<-merge(df,var2,by="seq",all=TRUE)
  #print(head(df))
}

df[is.na(df)] <- 0
write.table(df,paste0(opt$path,opt$otu),sep="\t",quote=FALSE)

#Get DNA seqs and make fasta file (sequences were already aligned by singleM)
dna.df<-df[,1]
dna.df<-as.data.frame(dna.df)
dna.df$'seq.name'<-row.names(dna.df)
colnames(dna.df)<-c("seq.text","seq.name")
fasta.df = data.frame(matrix("", ncol = 2, nrow = dim(df[1])))
fasta.df[c(1)]<-as.character(dna.df$seq.text)
fasta.df[c(2)]<-as.character(dna.df$seq.name)
colnames(fasta.df)<-c("seq.text","seq.name")

dat2fasta(fasta.df,outfile=paste0(opt$path,"otu.fasta"))
ex.fasta <- read.dna(paste0(path,"otu.fasta"),format="fasta")

#Cluster sequences in fasta file and keep one representatitive sequence per cluster
s.seed(9999)
otus<-otu(ex.fasta,k=5,threshold=as.numeric(opt$threshold),method="centroid")
clusterdf<-as.data.frame(cbind(otus,names(otus),dna.df))
clusterdf<-clusterdf[c(1,2,3)]
colnames(clusterdf)<-c("otu_num","cluster_rep","seq")
clustersub<-subset(clusterdf,endsWith(clusterdf$cluster_rep,suffix=c('*'))==TRUE)
clustersub<-clustersub[c(1,3)]
colnames(clustersub)<-c("otu_num","rep_seq")

#Aggregate counts by cluster
rownames(df)<-df$seq
df$seq<-NULL
df$otu_num<-otus
df_aggregate<-aggregate(.~otu_num,df,sum)
df_merged<-merge(df_aggregate,clustersub,by="otu_num")
write.table(df_merged,paste0(opt$path,opt$cluster),sep="\t",quote=FALSE)

