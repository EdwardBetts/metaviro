library(MASS)
library(caret)
library(reshape2)
library(ggplot2)
library(dplyr)
setwd("~/Documents/metaviro/experiments/140613_more_patterns/")

# we load the k-merized sequences as well as the estimated champions

input_sequences_f="long_kmers_m300_n4.fa_111.csv"
pattern_name=strsplit(input_sequences_f,"_|\\.")[[1]][6]
input_sequences=fread(input_sequences_f)
processed_k_mers=fread(paste("processed/",input_sequences_f,"_retrieval_stats.csv",sep=""))


s                                             

champions=unique(rbind(
  head(long_kmer_retrieval_d[!is.na(arch_f1_score)][order(arch_f1_score,decreasing=T)],n=30) %.% mutate("best"="A"),
  head(long_kmer_retrieval_d[!is.na(bact_f1_score)][order(bact_f1_score,decreasing=T)],n=30) %.% mutate("best"="B"),
  head(long_kmer_retrieval_d[!is.na(euk_f1_score)][order(euk_f1_score,decreasing=T)],n=30) %.% mutate("best"="E"),
  head(long_kmer_retrieval_d[!is.na(virus_f1_score)][order(virus_f1_score,decreasing=T)],n=30) %.% mutate("best"="V")
))


# plot in 2d space 

input_selected=input_sequences[kmer %in% champions$kmer]
count_m=dcast.data.table(input_selected,sequence_description+class+species~kmer,value.var='count',fill=0)
count_cols=colnames(count_m)[-(1:3)]

d <- dist(count_m[,count_cols,with=F]) # euclidean distances between the rows
fit <- cmdscale(d,eig=TRUE, k=4) # k is the number of dim
#fit <- isoMDS(d, k=2) # k is the number of dim
count_m$xcoord=fit$points[,1]
count_m$ycoord=fit$points[,2]
count_m$zcoord=fit$points[,3]
count_m$ucoord=fit$points[,4]

ggplot(count_m,aes(x=xcoord,y=ycoord,colour=class))+geom_point()
ggplot(count_m,aes(x=xcoord,y=zcoord,colour=class))+geom_point()
ggplot(count_m,aes(x=xcoord,y=ucoord,colour=class))+geom_point()
ggplot(count_m,aes(x=ycoord,y=zcoord,colour=class))+geom_point()
ggplot(count_m,aes(x=ycoord,y=ucoord,colour=class))+geom_point()

# Try tree reg 

set.seed(3456)
trainIndex <- createDataPartition(count_m$class, p = .8, list = FALSE, times = 1)

seqTrain <- count_m[ trainIndex[,1]][,c("class",count_cols),with=F]
seqTest  <- count_m[-trainIndex[,1],c("class",count_cols),with=F]

fitControl <- trainControl(    method = "repeatedcv", number = 3, repeats = 2)
rforestFit <- train(class ~ ., data = seqTrain,
                 method = "rf",
                 trControl = fitControl,
                 ## This last option is actually one
                 ## for gbm() that passes through
                 verbose = T)

gbmFit1
