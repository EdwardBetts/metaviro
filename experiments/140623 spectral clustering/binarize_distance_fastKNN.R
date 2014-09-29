library(gridExtra)
library(fields)
library(igraph)
library(corrplot)
library(ggplot2)
library(data.table)
library(reshape2)
library(plyr)
library(dplyr)

library(foreach)
library(doParallel)
n_cores=4
registerDoParallel(cores=n_cores)

# setwd("C:/Users/hayssam/Documents/GitHub/metaviro/experiments/140623 spectral clustering")
setwd("~/metaviro/experiments/140623 spectral clustering")
input_pattern="11111"
mixed_3mers = fread(paste('../140613_more_patterns/spaced_kmers/long_kmers_mVar_n30.fa_',input_pattern,".csv",sep=""))

mixed_3mers = dcast.data.table(mixed_3mers,pattern+GC+species+sequence_description~kmer,value.var='count',fill=0)

seq_attributes=limma::strsplit2(mixed_3mers$sequence_description,"_")
mixed_3mers$class=factor(seq_attributes[,2])
mixed_3mers$species=factor(seq_attributes[,1])
mixed_3mers$pattern=factor(mixed_3mers$pattern)

mixed_3mers$GC_f=cut(mixed_3mers$GC, quantile(mixed_3mers$GC,prob=seq(0,1,0.1)),include.lowest=T)

ggplot(mixed_3mers,aes(x=class,y=GC))+geom_boxplot()
ggplot(mixed_3mers,aes(x=class,y=GC))+geom_point()

setkey(mixed_3mers,"sequence_description")

mixed_3mers_label_cols=c('GC_f','class',"pattern",'GC','species','sequence_description')
mixed_3mers_count_cols=colnames(mixed_3mers)
mixed_3mers_count_cols=setdiff(mixed_3mers_count_cols,mixed_3mers_label_cols)


#with length normalization
rowSums(mixed_3mers[,mixed_3mers_count_cols,with=F])
mixed_3mers[16778,]

mixed_3mers_counts_scaled=scale(mixed_3mers_counts)


compute_closest= function(feat_matrix,selected_features,selected_obs,n_closer,selected_labels=NA){
	selected_obs=na.omit(selected_obs)
	dist_m=rdist(feat_matrix[selected_obs,selected_features,with=F],feat_matrix[,selected_features,with=F])
	ranked = apply(dist_m,MARGIN=1,function(r)rank(r,ties.method="min"))
	if(is.na(selected_labels)){
		ldply(1:length(selected_obs),function(i){	
			cat(i,"\n")
			selected_i=which(ranked[,i]<=n_closer); 
			data.table(src=selected_obs[i],tgt=selected_i,rank=ranked[selected_i,i],dist=dist_m[i,selected_i])
		})

	}else{
		ldply(1:length(selected_obs),function(i){	
			selected_i=which(ranked[,i]<=n_closer); 
			tgt_attributes=feat_matrix[selected_i,selected_labels,with=F];
			src_attributes=feat_matrix[selected_obs[i],selected_labels,with=F];
			setnames(src_attributes,paste("src",selected_labels,sep="_"));
			setnames(tgt_attributes,paste("tgt",selected_labels,sep="_"));
			cbind(src_attributes,data.table(src=selected_obs[i],tgt=selected_i,rank=ranked[selected_i,i],dist=dist_m[i,selected_i]),tgt_attributes)
		})
	}
}


input_size=nrow(mixed_3mers)
nas_to_add=((input_size %/% n_cores)+1)*n_cores - input_size
input_splits= matrix(c(1:input_size,rep(NA,nas_to_add)),byrow=T,ncol=n_cores)

res= data.table(foreach(i=iter(input_splits,by='col'),.combine='rbind',.inorder=F) %dopar% compute_closest(mixed_3mers,mixed_3mers_count_cols, i[,1],n_closer=10,selected_labels=c("class","species","sequence_description")))
save(res,file=paste(input_pattern,"closest",10,"dist_matrix.RData",sep="_"))

## Display the counts for the seq having a lot of close neighbors 
res[,.N,by=list(src,src_sequence_description,src_species)][order(N)] %.% tail()
res[,.N,by=list(tgt,tgt_sequence_description,tgt_species)][order(N)] %.% tail()

mixed_3mers[res[src==10096,tgt_sequence_description]]
mixed_3mers[res[src==10099,tgt_sequence_description]]




# compare distribution of intra-class vs inter-class

all_mean_dist=res[src_sequence_description!=tgt_sequence_description,mean(dist),by=list(src_class,tgt_class)]
all_mean_dist_m=dcast(all_mean_dist,src_class~tgt_class)
rownames(all_mean_dist_m)=all_mean_dist_m$src_class
all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
corrplot(all_mean_dist_m,is.corr=F)

ggplot(res,aes(x=tgt_class,colour=tgt_class,y=dist))+geom_point()+facet_wrap(~src_class)

all_mean_rank=res[src_sequence_description!=tgt_sequence_description,mean(rank),by=list(src_class,tgt_class)]
all_mean_rank_m=dcast(all_mean_rank,src_class~tgt_class)
rownames(all_mean_rank_m)=all_mean_rank_m$src_class
all_mean_rank_m=as.matrix(all_mean_rank_m[,-c(1)])
corrplot(all_mean_rank_m,is.corr=F)


#compare distributions by species
all_mean_dist=res[src_sequence_description!=tgt_sequence_description,mean(dist),by=list(src_species,tgt_species)]
all_mean_dist_m=dcast(all_mean_dist,src_species~tgt_species)
rownames(all_mean_dist_m)=all_mean_dist_m$src_species
all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
all_mean_dist_m[is.na(all_mean_dist_m)] <- max(all_mean_dist_m,na.rm=T)*2
corrplot(all_mean_dist_m[1:75,1:75],is.corr=F,order="hclust")

# using a graph to model 5 closest neighbors 
n_neighbors=2
pairwise_dist_g=graph.data.frame(res[rank<=n_neighbors,list(src_sequence_description,tgt_sequence_description)])
V(pairwise_dist_g)$class=as.character(mixed_3mers[V(pairwise_dist_g)$name,class][,class])
V(pairwise_dist_g)$species=as.character(mixed_3mers[V(pairwise_dist_g)$name,species][,species])
V(pairwise_dist_g)$GC=mixed_3mers[V(pairwise_dist_g)$name,GC][,GC]
V(pairwise_dist_g)$GC_f=as.character(mixed_3mers[V(pairwise_dist_g)$name,GC_f][,GC_f])

pairwise_dist_g_simple=simplify(pairwise_dist_g)
head(get.data.frame(pairwise_dist_g_simple,what='vertices'))

## degree distrib 

tail(sort(degree(pairwise_dist_g_simple)))
lines(density(degree(pairwise_dist_g_simple)))
write.graph(pairwise_dist_g_simple,format="gml",file=paste("dist_graph_mVar",input_pattern,"_",n_neighbors,".gml",sep=""))
#degree distrib ? 

