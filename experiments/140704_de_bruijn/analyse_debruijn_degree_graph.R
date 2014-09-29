library(ggplot2)
library(data.table)
library(corrplot)
library(reshape2)
library(doMC)

registerDoMC(cores=4)


input_counts=fread("sparse_matrix_5mers_10k_contigs.csv")

input_counts_m=dcast.data.table(input_counts,class+sequence_description~kmer,fill=0)



deg_by_class=input_counts[,mean(degree),by=list(class,kmer)]
deg_by_class_m=dcast.data.table(deg_by_class,kmer~ class)



ggplot(deg_by_class,aes(x=kmer,y=class,colour=V1))+geom_point()

deg_by_class[,ranked_by_deg:=rank(-1*V1),by=list(class)]
deg_by_class[ranked_by_deg<=5][order(kmer)]

dcast(deg_by_class[ranked_by_deg<=10],kmer~class,value.var='V1')
dcast(deg_by_class[ranked_by_deg<=20],kmer~class,value.var='V1')
dcast(deg_by_class[ranked_by_deg<=30],kmer~class,value.var='V1')
dcast(deg_by_class[ranked_by_deg<=100],kmer~class,value.var='V1')

dcast(deg_by_class[ranked_by_deg<=10],kmer~class,value.var='ranked_by_deg')
dcast(deg_by_class[ranked_by_deg<=100],kmer~class,value.var='ranked_by_deg')


input_counts_m[,list(class,CGTCG)][,mean(CGTCG),by=class]
input_counts_m[,list(class,GAAGA)][,mean(GAAGA),by=class]
input_counts_m[,list(class,AAGAA)][,mean(AAGAA),by=class]
input_counts_m[,list(class,AAAGA)][,mean(AAAGA),by=class]

## Frequencies of each k-mer being used as a branch point 


r=input_counts_m[,list(class,CGTCG)]
by(r,r$class,function(s)table(s$CGTCG>0))
by(r,r$class,function(s)table(s$CGTCG>3))

r=input_counts_m[,list(class,TGATG)]
by(r,r$class,function(s)table(s$TGATG>0))
by(r,r$class,function(s)table(s$TGATG>3))


#
r=input_counts_m[,list(class,GCGGC)]
xtabs(~class+(GCGGC>6),r)


# 
r=input_counts_m[,list(class,GAAGA)]
xtabs(~class+(GAAGA>3),r)

r=input_counts_m[,list(class,GAAAA)]
xtabs(~class+(GAAAA>3),r)


# bacterian break point 
r=input_counts_m[,list(class,CGTCG)]
xtabs(~class+(CGTCG>6),r)


# we subset with the high degree k-mers and try a tree classifier

high_degree_kmers=dcast(deg_by_class[ranked_by_deg<=60],kmer~class,value.var='V1')$kmer
selected_kmers=input_counts_m[,c("class",high_degree_kmers),with=F]
selected_kmers$class=factor(selected_kmers$class)

# Remove zero/1 everywhere contig ? 

selected_kmers=selected_kmers[rowSums(selected_kmers[,high_degree_kmers,with=F])>=10]

library(caret)
library(party)
inTrain=createDataPartition(selected_kmers$class,p=0.5,list=F)[,1]

trainFeat=selected_kmers[inTrain,high_degree_kmers,with=F]
trainClass=selected_kmers[inTrain,class]
testFeat=selected_kmers[-inTrain,high_degree_kmers,with=F]
testClass=selected_kmers[-inTrain,class]

prop.table(table(trainClass))
prop.table(table(testClass))

# any correlations between them ? 
c orrplot(cor(trainFeat))
svmFit = train(trainFeat,trainClass,method="svmRadial",metric="Kappa",verbose=T,trControl=trainControl(method="CV",number=4),tuneGrid=expand.grid(sigma=c(0.006,1/90),C=1))

confusionMatrix(testClass,predict(svmFit,testFeat))


# Rforest 
rforestFit=train(trainFeat,trainClass,method="rf",metric="Kappa",verbose=T,trControl=trainControl(method="CV",number=4))
confusionMatrix(testClass,predict(rforestFit,testFeat))