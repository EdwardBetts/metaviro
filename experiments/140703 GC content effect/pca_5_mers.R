library(ggplot2)
require(MASS)
library(caret)
library(GGally)
library(dplyr)
library(plyr)
library(data.table)
library(reshape2)

setwd("/Users/hayssam/Documents/metaviro/experiments/140703 GC content effect/")
all_kmers=fread("../140613_more_patterns/spaced_kmers/long_kmers_mVar_n30.fa_11111.csv")

all_kmers_counts = dcast.data.table(all_kmers,sequence_description+GC+sequence_length+species~kmer,value.var="count",fill=0)
all_kmers_counts$class=unlist(lapply(strsplit(all_kmers_counts$sequence_description,"_"),"[",2))
count_columns=colnames(all_kmers_counts)
count_columns = setdiff(count_columns,c("sequence_description","GC","sequence_length","species","class"))

# We first normalize by contig lenght 
all_kmers_counts_normalized=as.data.table(apply(all_kmers_counts[,count_columns,with=F],MARGIN=2,FUN=function(l) l/all_kmers_counts$sequence_length))
all_kmers_counts_normalized=cbind(all_kmers_counts_normalized,all_kmers_counts[,list(class)])


# We use the caret tranformations to do a PCA 
trans = preProcess(all_kmers_counts_normalized[,count_columns,with=F], 
                   method=c("BoxCox", "center", 
                            "scale", "pca"))
PC = predict(trans, all_kmers_counts_normalized[,count_columns,with=F])
dim(PC)
library(ggbiplot)


all_kmers_counts=cbind(all_kmers_counts,PC[,1:100])
pca_coord_cols=colnames(PC)[1:100]
all_kmers_counts=all_kmers_counts[,unique(colnames(all_kmers_counts)),with=F]


ggplot(all_kmers_counts,aes(x=PC1,y=PC2,colour=class))+geom_point()
ggplot(all_kmers_counts,aes(x=PC1,y=PC2,colour=GC,size=sequence_length))+geom_point()+ggtitle("First 2 PC on unscaled 4 mers data")
ggplot(all_kmers_counts,aes(x=PC2,y=PC3,colour=class))+geom_point()
ggplot(all_kmers_counts,aes(x=PC2,y=PC3,colour=class))+geom_point()+facet_wrap(~class)

ggsave(file="First 2 PC on boxcox, center, scaled and PCA'ed 5 mers data.pdf")



# let's compute the kMeans of each class 

kmeans_by_classes=rbind(
	data.table(kmeans(all_kmers_counts[class=="bact",pca_coord_cols,with=F],centers=20)$centers,class="bact"),
	data.table(kmeans(all_kmers_counts[class=="euk",pca_coord_cols,with=F],centers=20)$centers,class="euk"),
	data.table(kmeans(all_kmers_counts[class=="archea",pca_coord_cols,with=F],centers=20)$centers,class="archea"),
	data.table(kmeans(all_kmers_counts[class=="viruses",pca_coord_cols,with=F],centers=20)$centers,class="viruses"))



ggplot(all_kmers_counts,aes(x=PC1,y=PC2,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC1,y=PC2,fill=class),colour="black")
ggplot(all_kmers_counts,aes(x=PC2,y=PC3,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC2,y=PC3,fill=class),colour="black")
ggplot(all_kmers_counts,aes(x=PC2,y=PC4,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC2,y=PC4,fill=class),colour="black")
ggplot(all_kmers_counts,aes(x=PC3,y=PC4,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3,y=PC4,fill=class),colour="black")
ggplot(all_kmers_counts,aes(x=PC3,y=PC5,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3,y=PC5,fill=class),colour="black")
ggplot(all_kmers_counts,aes(x=PC4,y=PC5,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC4,y=PC5,fill=class),colour="black")


## Trying a kNN 15 with caret 
library(doMC)
registerDoMC(4)
inTrain <- createDataPartition(all_kmers_counts$class, p = 0.8, list = FALSE)

trainExpr <- as.matrix(all_kmers_counts[inTrain[,1],pca_coord_cols,with=F])
testExpr <- as.matrix(all_kmers_counts[-inTrain[,1],pca_coord_cols,with=F])

trainClass <- factor(all_kmers_counts[inTrain[,1],class])
testClass <- factor(all_kmers_counts[-inTrain[,1],class])

cv_opts = trainControl(method="cv", number=2)
knn_opts = data.frame(.k=c(11,13,15,17,19)) #odd to avoid ties

knn_performances <- train(x = trainExpr, y = trainClass,method = "knn",trControl=cv_opts,tuneGrid=knn_opts)

predicted_classes=predict(knn_performances,testExpr)
confusionMatrix(predicted_classes,testClass)

# Vs an SVM-PCA
svm_opts = expand.grid(.sigma=c(1/80,1/90,1/100,1/110),.C=c(0.1,1,10,20)) #odd to avoid ties
svm_performances <- train(x = trainExpr, y = trainClass,method = "svmRadial",trControl=cv_opts,tuneGrid=svm_opts,verbose=T)
predicted_classes=predict(svm_performances,testExpr)
confusionMatrix(predicted_classes,testClass)
save(svm_performances,trans,all_kmers_counts,file="svm-rbf-pca-90.Rdat")


# if we drop the first PC 
cols_to_keep=pca_coord_cols[-c(1)]
knn_performances_noPC1 <- train(x = trainExpr[,cols_to_keep], y = trainClass,method = "knn",trControl=cv_opts,tuneGrid=knn_opts)
predicted_classes_no_PC1=predict(knn_performances_noPC1,testExpr[,cols_to_keep])
confusionMatrix(predicted_classes_no_PC1,testClass)


# Vs the full matrix 
trainExpr <- as.matrix(all_kmers_counts[inTrain[,1],count_columns,with=F])
testExpr <- as.matrix(all_kmers_counts[-inTrain[,1],count_columns,with=F])
knn_performances_Full <- train(x = trainExpr, y = trainClass,method = "knn",trControl=cv_opts,tuneGrid=knn_opts)




# Trying a nnet 

nnet_perf <- train(x = trainExpr[,pca_coord_cols], y = trainClass,method = "nnet",trControl=cv_opts)
predicted_classes_nnet=predict(nnet_perf,testExpr[,pca_coord_cols])
confusionMatrix(predicted_classes_nnet,testClass)
