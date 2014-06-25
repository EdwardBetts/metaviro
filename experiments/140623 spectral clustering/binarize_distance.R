library(igraph)
library(corrplot)
library(ggplot2)
library(data.table)
library(reshape2)
setwd("C:/Users/hayssam/Documents/GitHub/metaviro/experiments/140623 spectral clustering")

mixed_3mers = fread('../140613_more_patterns/spaced_kmers/long_kmers_m300_n4.fa_111.csv')
seq_attributes=limma::strsplit2(mixed_3mers$sequence_description,"_")
mixed_3mers$class=factor(seq_attributes[,2])
mixed_3mers$species=factor(seq_attributes[,1])
mixed_3mers$pattern=factor(mixed_3mers$pattern)
mixed_3mers$kmer=factor(mixed_3mers$kmer)

mixed_3mers = dcast.data.table(mixed_3mers,GC+class+species+sequence_description~kmer,value.var='count',fill=0)

mixed_3mers$GC_f=cut(mixed_3mers$GC, quantile(mixed_3mers$GC,prob=seq(0,1,0.1)),include.lowest=T)

ggplot(mixed_3mers,aes(x=class,y=GC))+geom_boxplot()
ggplot(mixed_3mers,aes(x=class,y=GC))+geom_point()

setkey(mixed_3mers,"sequence_description")

mixed_3mers_count_cols=colnames(mixed_3mers)[-c(1:4)]
mixed_3mers_count_cols=mixed_3mers_count_cols[mixed_3mers_count_cols!='GC_f']
mixed_3mers_counts=as.matrix(mixed_3mers[,mixed_3mers_count_cols,with=F])
head(mixed_3mers_counts)

pairwise_dist=dist(mixed_3mers_counts)
dist_m=as.matrix(pairwise_dist)
colnames(dist_m)=mixed_3mers[,sequence_description]
rownames(dist_m)=mixed_3mers[,sequence_description]

pairwise_dist_df=data.table(melt(dist_m))
setnames(pairwise_dist_df,c("src","tgt","dist"))
pairwise_dist_df_ranked=pairwise_dist_df[,ranks:=rank(dist),by=src]
pairwise_dist_df_ranked[ranks<=5][order(src)]



ggplot(melt(mixed_3mers[names(head(sort(dist_m[100,])))]),aes(x=variable,y=value,group=sequence_description,colour=sequence_description))+geom_point()+geom_line()


# compare distribution of intra-class vs inter-class

colnames(dist_m)=mixed_3mers[,class]
rownames(dist_m)=mixed_3mers[,class]
dist_classes=data.table(melt(dist_m))

all_mean_dist=dist_classes[,mean(value,na.rm=T),by=list(Var1,Var2)]
all_mean_dist_m=dcast(all_mean_dist,Var1~Var2)
rownames(all_mean_dist_m)=all_mean_dist_m$Var1
all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
corrplot(all_mean_dist_m,is.corr=F)


#compare distributions by species
colnames(dist_m)=mixed_3mers[,species]
rownames(dist_m)=mixed_3mers[,species]
dist_species=data.table(melt(dist_m))

all_mean_dist_species=dist_species[,mean(value,na.rm=T),by=list(Var1,Var2)]

all_mean_dist_m=dcast(all_mean_dist_species,Var1~Var2)
rownames(all_mean_dist_m)=all_mean_dist_m$Var1
all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
corrplot(all_mean_dist_m[1:50,1:50],is.corr=F,order="hclust")


mixed_3mers[species=="10937870"]
mixed_3mers[species=="336279"]

# compare by GC class


colnames(dist_m)=mixed_3mers[,GC_f]
rownames(dist_m)=mixed_3mers[,GC_f]
dist_GC=data.table(melt(dist_m))

all_mean_dist_GC=dist_GC[,mean(value,na.rm=T),by=list(Var1,Var2)]

all_mean_dist_m=dcast(all_mean_dist_GC,Var1~Var2)
rownames(all_mean_dist_m)=all_mean_dist_m$Var1
all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
corrplot(all_mean_dist_m,is.corr=F)



# using a graph to model 5 closest neighbors 

pairwise_dist_g=graph.data.frame(pairwise_dist_df_ranked[ranks<=5])
V(pairwise_dist_g)$class=as.character(mixed_3mers[V(pairwise_dist_g)$name,class][,class])
V(pairwise_dist_g)$species=as.character(mixed_3mers[V(pairwise_dist_g)$name,species][,species])
V(pairwise_dist_g)$GC=mixed_3mers[V(pairwise_dist_g)$name,GC][,GC]
V(pairwise_dist_g)$GC_f=as.character(mixed_3mers[V(pairwise_dist_g)$name,GC_f][,GC_f])
head(get.data.frame(pairwise_dist_g,what='vertices'))

#degree distrib ? 

