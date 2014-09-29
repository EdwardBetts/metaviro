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

pc_kmers=princomp(all_kmers_counts[,count_columns,with=F])

coord_pca=predict(pc_kmers)[,1:3]

all_kmers_counts=cbind(all_kmers_counts,coord_pca)

ggplot(all_kmers_counts,aes(x=Comp.1,y=Comp.2,colour=class))+geom_point()
ggplot(all_kmers_counts,aes(x=Comp.1,y=Comp.2,colour=GC,size=sequence_length))+geom_point()+ggtitle("First 2 PC on unscaled 4 mers data")
ggsave(file="First 2 PC on unscaled 4 mers data.pdf")


## Same with scaling

pc_kmers=prcomp(all_kmers_counts[,count_columns,with=F],scale=T)
coord_pca=predict(pc_kmers)[,1:3]
colnames(coord_pca) <- paste(colnames(coord_pca),"_scaled",sep="")
all_kmers_counts=cbind(all_kmers_counts,coord_pca)

ggplot(all_kmers_counts,aes(x=PC1_scaled,y=PC2_scaled,colour=class))+geom_point()
ggplot(all_kmers_counts,aes(x=PC1_scaled,y=PC2_scaled,colour=GC,size=sequence_length))+geom_point()+ggtitle("First 2 PC on unscaled 4 mers data")


# Doesn't help 

## Same with length normalization 

all_kmers_counts_normalized=apply(all_kmers_counts[,count_columns,with=F],MARGIN=2,FUN=function(l) l/all_kmers_counts$sequence_length)
pc_kmers=prcomp(all_kmers_counts_normalized,scale=T)
coord_pca=predict(pc_kmers)[,1:10]
colnames(coord_pca) <- paste(colnames(coord_pca),"_scaled",sep="")
all_kmers_counts_normalized=data.table(cbind(all_kmers_counts_normalized,all_kmers_counts[,list(class,GC,sequence_description,sequence_length)],coord_pca))


ggplot(all_kmers_counts_normalized,aes(x=PC1_scaled,y=PC2_scaled,colour=class))+geom_point()
ggplot(all_kmers_counts_normalized,aes(x=PC1_scaled,y=PC2_scaled,colour=GC,size=sequence_length))+geom_point()+ggtitle("First 2 PC on scaled 4 mers data")
ggplot(all_kmers_counts_normalized,aes(x=PC1_scaled,y=PC2_scaled,colour=class,size=GC))+geom_point()+ggtitle("First 2 PC on scaled 4 mers data")
ggplot(all_kmers_counts_normalized,aes(x=PC2_scaled,y=PC3_scaled,colour=class,size=GC))+geom_point(alpha=0.3)+ggtitle("First 2 PC on scaled 4 mers data")+geom_density2d()



# let's compute the kMeans of each class 

kmeans_by_classes=rbind(
	data.table(kmeans(all_kmers_counts_normalized[class=="bact",colnames(coord_pca),with=F],centers=20)$centers,class="bact"),
	data.table(kmeans(all_kmers_counts_normalized[class=="euk",colnames(coord_pca),with=F],centers=20)$centers,class="euk"),
	data.table(kmeans(all_kmers_counts_normalized[class=="archea",colnames(coord_pca),with=F],centers=20)$centers,class="archea"),
	data.table(kmeans(all_kmers_counts_normalized[class=="viruses",colnames(coord_pca),with=F],centers=20)$centers,class="viruses"))



ggplot(all_kmers_counts_normalized,aes(x=PC1_scaled,y=PC2_scaled,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC1_scaled,y=PC2_scaled,fill=class),colour="black")

ggplot(all_kmers_counts_normalized,aes(x=PC3_scaled,y=PC2_scaled,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3_scaled,y=PC2_scaled,fill=class,colour=class))
ggplot(all_kmers_counts_normalized,aes(x=PC4_scaled,y=PC2_scaled,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3_scaled,y=PC2_scaled,fill=class,colour=class))
ggplot(all_kmers_counts_normalized,aes(x=PC5_scaled,y=PC3_scaled,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3_scaled,y=PC2_scaled,fill=class,colour=class))
ggplot(all_kmers_counts_normalized,aes(x=PC5_scaled,y=PC4_scaled,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3_scaled,y=PC2_scaled,fill=class,colour=class))
ggplot(all_kmers_counts_normalized,aes(x=PC6_scaled,y=PC5_scaled,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3_scaled,y=PC2_scaled,fill=class,colour=class))


## Trying PCA on a GC-bounded range of species 

med_GC=all_kmers_counts[GC>=20 & GC<=60]
med_GC_normalized=apply(med_GC[,count_columns,with=F],MARGIN=2,FUN=function(l) l/med_GC$sequence_length)
pc_kmers=prcomp(med_GC_normalized,scale=T)
coord_pca=predict(pc_kmers)[,1:20]
colnames(coord_pca) <- paste(colnames(coord_pca),"_scaled",sep="")
med_GC_normalized=data.table(cbind(med_GC_normalized,med_GC[,list(class,GC,sequence_description,sequence_length)],coord_pca))

kmeans_by_classes=rbind(
	data.table(kmeans(med_GC_normalized[class=="bact",colnames(coord_pca),with=F],centers=20)$centers,class="bact"),
	data.table(kmeans(med_GC_normalized[class=="euk",colnames(coord_pca),with=F],centers=20)$centers,class="euk"),
	data.table(kmeans(med_GC_normalized[class=="archea",colnames(coord_pca),with=F],centers=20)$centers,class="archea"),
	data.table(kmeans(med_GC_normalized[class=="viruses",colnames(coord_pca),with=F],centers=20)$centers,class="viruses"))


ggplot(med_GC_normalized,aes(x=PC1_scaled,y=PC2_scaled,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC1_scaled,y=PC2_scaled,fill=class,colour=class))

ggplot(med_GC_normalized,aes(x=PC3,y=PC2,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3,y=PC2,fill=class,colour=class))
ggplot(med_GC_normalized,aes(x=PC4,y=PC2,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3,y=PC2,fill=class,colour=class))
ggplot(med_GC_normalized,aes(x=PC4,y=PC3,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=10,data=kmeans_by_classes,aes(x=PC3,y=PC2,fill=class,colour=class))



## Trying a kNN 15 with caret 

inTrain <- createDataPartition(med_GC_normalized$class, p = 0.8, list = FALSE)

trainExpr <- as.matrix(med_GC_normalized[inTrain[,1],colnames(coord_pca),with=F])
testExpr <- as.matrix(med_GC_normalized[-inTrain[,1],colnames(coord_pca),with=F])

trainClass <- factor(med_GC_normalized[inTrain[,1],class])
testClass <- factor(med_GC_normalized[-inTrain[,1],class])

cv_opts = trainControl(method="cv", number=2)
knn_opts = data.frame(.k=c(11,13,15,17,19)) #odd to avoid ties

knn_performances <- train(x = trainExpr, y = trainClass,method = "knn",trControl=cv_opts,tuneGrid=knn_opts)

predicted_classes=predict(knn_performances,testExpr)
confusionMatrix(predicted_classes,testClass)


# if we drop the first PC 
cols_to_keep=colnames(coord_pca)[-c(1)]
trainExpr <- as.matrix(med_GC_normalized[inTrain[,1],cols_to_keep,with=F])
testExpr <- as.matrix(med_GC_normalized[-inTrain[,1],cols_to_keep,with=F])

trainClass <- factor(med_GC_normalized[inTrain[,1],class])
testClass <- factor(med_GC_normalized[-inTrain[,1],class])

cv_opts = trainControl(method="cv", number=2)
knn_opts = data.frame(.k=c(11,13,15,17,19,21)) #odd to avoid ties

knn_performances_noPC1 <- train(x = trainExpr, y = trainClass,method = "knn",trControl=cv_opts,tuneGrid=knn_opts)

predicted_classes=predict(knn_performances_noPC1,testExpr)
confusionMatrix(predicted_classes,testClass)




# Versus knn on the full matrix 
# performances do drop

cols_to_keep=count_columns
trainExpr <- as.matrix(med_GC_normalized[inTrain[,1],cols_to_keep,with=F])
testExpr <- as.matrix(med_GC_normalized[-inTrain[,1],cols_to_keep,with=F])

trainClass <- factor(med_GC_normalized[inTrain[,1],class])
testClass <- factor(med_GC_normalized[-inTrain[,1],class])

cv_opts = trainControl(method="cv", number=2)
knn_opts = data.frame(.k=c(15)) #odd to avoid ties

knn_performances_full <- train(x = trainExpr, y = trainClass,method = "knn",trControl=cv_opts,tuneGrid=knn_opts)

predicted_classes_full=predict(knn_performances_full,testExpr)
confusionMatrix(predicted_classes_full,testClass)



# Trying with more dimensions 





coord_pca_more=predict(pc_kmers)[,21:32]
colnames(coord_pca_more) <- paste(colnames(coord_pca_more),"_scaled",sep="")
med_GC_normalized=cbind(med_GC_normalized,coord_pca_more)
coord_pca_cols = grep("PC",colnames(med_GC_normalized),value=T)

trainExpr <- as.matrix(med_GC_normalized[inTrain[,1],coord_pca_cols,with=F])
testExpr <- as.matrix(med_GC_normalized[-inTrain[,1],coord_pca_cols,with=F])

trainClass <- factor(med_GC_normalized[inTrain[,1],class])
testClass <- factor(med_GC_normalized[-inTrain[,1],class])

cv_opts = trainControl(method="cv", number=2)
knn_opts = data.frame(.k=c(11,13,15,17,19)) #odd to avoid ties

knn_performances <- train(x = trainExpr, y = trainClass,method = "knn",trControl=cv_opts,tuneGrid=knn_opts)

predicted_classes=predict(knn_performances,testExpr)
confusionMatrix(predicted_classes,testClass)



# We do a PCA that is depdendent on the GC content 

all_kmers_counts$GC_f=cut(all_kmers_counts$GC,quantile(all_kmers_counts$GC,probs=c(seq(0,1,0.3),1)),include.lowest=T)
all_kmers_counts_normalized=as.data.table(apply(all_kmers_counts[,count_columns,with=F],MARGIN=2,FUN=function(l) l/all_kmers_counts$sequence_length))
all_kmers_counts_normalized=cbind(all_kmers_counts_normalized,all_kmers_counts[,list(class,GC_f,GC,sequence_length,sequence_description)])
summary(all_kmers_counts_normalized$GC_f)

gc_wise_pca=data.table(ldply(levels(all_kmers_counts_normalized$GC_f),function(GC_level){
	subset=all_kmers_counts_normalized[GC_f==GC_level]
	pca=prcomp(subset[,count_columns,with=F],scale=T)
	coord_pca=predict(pca)[,1:32]
	return(cbind(subset,coord_pca))
}))

pca_cols=grep("PC",colnames(gc_wise_pca),value=T)

kmeans_by_classes=data.table(ldply(levels(gc_wise_pca$GC_f),function(GC_level){
	these_kmeans=rbind(
		data.table(kmeans(gc_wise_pca[class=="bact",pca_cols,with=F],centers=10)$centers,class="bact"),
		data.table(kmeans(gc_wise_pca[class=="euk",pca_cols,with=F],centers=10)$centers,class="euk"),
		data.table(kmeans(gc_wise_pca[class=="archea",pca_cols,with=F],centers=10)$centers,class="archea"),
		data.table(kmeans(gc_wise_pca[class=="viruses",pca_cols,with=F],centers=10)$centers,class="viruses"))
	these_kmeans$GC_f=GC_level
	return(these_kmeans)
}))


ggplot(gc_wise_pca,aes(x=PC1,y=PC2,colour=class))+geom_point()
ggplot(gc_wise_pca,aes(x=PC1,y=PC2,colour=class))+geom_point(alpha=0.3)+facet_wrap(~GC_f)
ggplot(gc_wise_pca,aes(x=PC1,y=PC2,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=3,data=kmeans_by_classes,aes(x=PC1,y=PC2,fill=class),colour="black")+facet_wrap(~GC_f,scale="free")

ggplot(gc_wise_pca,aes(x=PC1,y=PC2,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=3,data=kmeans_by_classes,aes(x=PC1,y=PC2,fill=class),colour="black")+facet_grid(GC_f~class)
ggplot(gc_wise_pca,aes(x=PC3,y=PC2,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=3,data=kmeans_by_classes,aes(x=PC3,y=PC2,fill=class),colour="black")+facet_grid(GC_f~class)
ggplot(gc_wise_pca,aes(x=PC11,y=PC22,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=3,data=kmeans_by_classes,aes(x=PC11,y=PC22,fill=class),colour="black")+facet_grid(GC_f~class)

ggplot(gc_wise_pca,aes(x=PC3,y=PC2,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=3,data=kmeans_by_classes,aes(x=PC3,y=PC2,fill=class),colour="black")+facet_wrap(~GC_f,scale="free")

ggplot(gc_wise_pca,aes(x=PC22,y=PC11,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=3,data=kmeans_by_classes,aes(x=PC22,y=PC11,fill=class),colour="black")+facet_wrap(~GC_f,scale="free")

ggplot(gc_wise_pca,aes(x=PC22,y=PC11,colour=class))+geom_point(alpha=0.3)+geom_point(shape=24,size=3,data=kmeans_by_classes,aes(x=PC22,y=PC11,fill=class),colour="black")+facet_wrap(~GC_f,scale="free")


# We try a kNN on the GC wise PCA 

inTrain <- createDataPartition(gc_wise_pca$GC_f, p = 0.8, list = FALSE)

trainExpr <- as.matrix(gc_wise_pca[inTrain[,1],pca_cols,with=F])
testExpr <- as.matrix(gc_wise_pca[-inTrain[,1],pca_cols,with=F])

trainClass <- factor(gc_wise_pca[inTrain[,1],class])
testClass <- factor(gc_wise_pca[-inTrain[,1],class])

cv_opts = trainControl(method="cv", number=2)
knn_opts = data.frame(.k=c(11,15,19)) #odd to avoid ties

knn_performances <- train(x = trainExpr, y = trainClass,method = "knn",trControl=cv_opts,tuneGrid=knn_opts)

predicted_classes=predict(knn_performances,testExpr)
confusionMatrix(predicted_classes,testClass)


# Very bad performances with 10 cuts and 32 dimensions





## Trying an LDA 
all_kmers_counts$GC_f=cut(all_kmers_counts$GC,quantile(all_kmers_counts$GC,probs=c(seq(0,1,0.3),1)),include.lowest=T)
all_kmers_counts_normalized=as.data.table(apply(all_kmers_counts[,count_columns,with=F],MARGIN=2,FUN=function(l) l/all_kmers_counts$sequence_length))
all_kmers_counts_normalized=cbind(all_kmers_counts_normalized,all_kmers_counts[,list(class)])

inTrain <- createDataPartition(all_kmers_counts_normalized$class, p = 0.5, list = FALSE)
trainData=as.data.frame(all_kmers_counts_normalized[inTrain[,1]])
testData=as.data.frame(all_kmers_counts_normalized[-inTrain[,1]])

lda_model=lda(formula=class~.,data=trainData)
lda_coord=predict(lda_model)
trainData=cbind(trainData,lda_coord$x)


ggplot(trainData,aes(x=LD1,y=LD2,colour=class))+geom_point(alpha=0.3)
ggplot(trainData,aes(x=LD3,y=LD2,colour=class))+geom_point(alpha=0.3)
ggpairs(data=trainData,columns=258:260,colour="class")


plda=predict(lda_model,testData)
confusionMatrix(lda_coord$class,trainData$class)
confusionMatrix(plda$class,testData$class)