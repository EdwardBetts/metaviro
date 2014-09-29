library(caret)
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
input_pattern="1111"
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
normalized_counts=mixed_3mers[,mixed_3mers_count_cols,with=F] / rowSums(mixed_3mers[,mixed_3mers_count_cols,with=F])

normalized_counts=cbind(normalized_counts,mixed_3mers[,mixed_3mers_label_cols,with=F])

mixed_3mers_counts_scaled=scale(mixed_3mers[,mixed_3mers_count_cols,with=F])
mixed_3mers_counts_scaled=cbind(mixed_3mers_counts_scaled,mixed_3mers[,mixed_3mers_label_cols,with=F])


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


n_closer_neighbors=50
input_size=nrow(mixed_3mers)
nas_to_add=((input_size %/% n_cores)+1)*n_cores - input_size
input_splits= matrix(c(1:input_size,rep(NA,nas_to_add)),byrow=T,ncol=n_cores)

res= data.table(foreach(i=iter(input_splits,by='col'),.combine='rbind',.inorder=F) %dopar% compute_closest(mixed_3mers,mixed_3mers_count_cols, i[,1],n_closer=n_closer_neighbors,selected_labels=c("class","species","sequence_description")))
res$metric=factor('euclidean',levels=c('euclidean','length_norm','scaled'))
res$pattern=input_pattern
save(res,file=paste(input_pattern,"closest",n_closer_neighbors,"dist_matrix.RData",sep="_"))

res_length_norm=data.table(foreach(i=iter(input_splits,by='col'),.combine='rbind',.inorder=F) %dopar% compute_closest(normalized_counts,mixed_3mers_count_cols, i[,1],n_closer=n_closer_neighbors,selected_labels=c("class","species","sequence_description")))
res_length_norm$metric=factor('length_norm',levels=c('euclidean','length_norm','scaled'))
res_length_norm$pattern=input_pattern
save(res_length_norm,file=paste(input_pattern,"closest",n_closer_neighbors,"length_norm_dist_matrix.RData",sep="_"))

res_scaled=data.table(foreach(i=iter(input_splits,by='col'),.combine='rbind',.inorder=F) %dopar% compute_closest(mixed_3mers_counts_scaled,mixed_3mers_count_cols, i[,1],n_closer=n_closer_neighbors,selected_labels=c("class","species","sequence_description")))
res_scaled$metric=factor('scaled',levels=c('euclidean','length_norm','scaled'))
res_scaled$pattern=input_pattern

save(res_scaled,file=paste(input_pattern,"closest",n_closer_neighbors,"scaled_dist_matrix.RData",sep="_"))




## Display the counts for the seq having a lot of close neighbors 
res[,.N,by=list(src,src_sequence_description,src_species)][order(N)] %.% tail()
res_length_norm[,.N,by=list(src,src_sequence_description,src_species)][order(N)] %.% tail()

res[,.N,by=list(tgt,tgt_sequence_description,tgt_species)][order(N)] %.% tail()
res_length_norm[,.N,by=list(tgt,tgt_sequence_description,tgt_species)][order(N)] %.% tail()

mixed_3mers[res[src==10096,tgt_sequence_description]]
mixed_3mers[res[tgt==10857,src_sequence_description]]
mixed_3mers[res[src==n_neighbors,tgt_sequence_description]]




# compare distribution of intra-class vs inter-class
old_par=par(mfrow=c(3,1))
all_mean_dist=res[src_sequence_description!=tgt_sequence_description,mean(dist),by=list(src_class,tgt_class)]
all_mean_dist_m=dcast(all_mean_dist,src_class~tgt_class)
rownames(all_mean_dist_m)=all_mean_dist_m$src_class
all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
corrplot(all_mean_dist_m,is.corr=F)

all_mean_dist=res_length_norm[src_sequence_description!=tgt_sequence_description,mean(dist),by=list(src_class,tgt_class)]
all_mean_dist_m=dcast(all_mean_dist,src_class~tgt_class)
rownames(all_mean_dist_m)=all_mean_dist_m$src_class
all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
corrplot(all_mean_dist_m,is.corr=F)

all_mean_dist=res_scaled[src_sequence_description!=tgt_sequence_description,mean(dist),by=list(src_class,tgt_class)]
all_mean_dist_m=dcast(all_mean_dist,src_class~tgt_class)
rownames(all_mean_dist_m)=all_mean_dist_m$src_class
all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
corrplot(all_mean_dist_m,is.corr=F)

par(old_par)

sample_rows = function(dt,size){dt[sample(1:nrow(dt),size),]}
ggplot(sample_rows(rbind(res,res_length_norm,res_scaled),size=50000),aes(x=tgt_class,shape=tgt_class,y=dist,size=rank<=5))+geom_point()+facet_grid(metric~src_class,scale='free')


# compare distribution of intra-class vs inter-class ranks
{
	old_par=par(mfrow=c(2,2))
	all_mean_dist=res[src_sequence_description!=tgt_sequence_description,mean(rank),by=list(src_class,tgt_class)]
	all_mean_dist_m=dcast(all_mean_dist,src_class~tgt_class)
	rownames(all_mean_dist_m)=all_mean_dist_m$src_class
	all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
	corrplot(all_mean_dist_m,is.corr=F,main=paste('unscaled',input_pattern))

	all_mean_dist=res_length_norm[src_sequence_description!=tgt_sequence_description,mean(rank),by=list(src_class,tgt_class)]
	all_mean_dist_m=dcast(all_mean_dist,src_class~tgt_class)
	rownames(all_mean_dist_m)=all_mean_dist_m$src_class
	all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
	corrplot(all_mean_dist_m,is.corr=F,main=paste('length_norm',input_pattern))

	all_mean_dist=res_scaled[src_sequence_description!=tgt_sequence_description,mean(rank),by=list(src_class,tgt_class)]
	all_mean_dist_m=dcast(all_mean_dist,src_class~tgt_class)
	rownames(all_mean_dist_m)=all_mean_dist_m$src_class
	all_mean_dist_m=as.matrix(all_mean_dist_m[,-c(1)])
	corrplot(all_mean_dist_m,is.corr=F,main=paste('scaled',input_pattern))

	par(old_par)
}

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





# Determine % of correct classification as a function of the number of neighbors



get_statistics = function(mat,training_seq,label,n_neighbors){
	unscaled_kNN=sample_rows(mat,size=120000)[rank<=n_neighbors,names(sort(table(tgt_class))[4]),by=list(src,src_sequence_description,src_class)]
	these_stats=as.data.table(t(confusionMatrix(unscaled_kNN$src_class,unscaled_kNN$V1)$overall))
	these_stats$n_neighbors=n_neighbors
	these_stats$label=label
	return(these_stats)
}

stat_by_neigh=rbind(
	ldply(3:30,function(n)get_statistics(res,"unscaled",n),.parallel=T,.progress="tk"),
	ldply(3:30,function(n)get_statistics(res,"unscaled",n),.parallel=T,.progress="tk"),
	ldply(3:30,function(n)get_statistics(res,"unscaled",n),.parallel=T,.progress="tk"),
	ldply(3:30,function(n)get_statistics(res,"unscaled",n),.parallel=T,.progress="tk"))

ggplot(stat_by_neigh,aes(x=n_neighbors,fill=label,y=Kappa))+geom_point()+geom_smooth()

