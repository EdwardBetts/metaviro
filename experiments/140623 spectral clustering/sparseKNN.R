{ #setup
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
	setwd("~/metaviro/experiments/140623 spectral clustering")
}

compute_closest_train_test= function(feat_matrix,train_indices,test_indices,selected_features,n_closer,selected_labels=NULL){
	test_indices=na.omit(test_indices)
	train_indices=na.omit(train_indices)
	dist_m=rdist(feat_matrix[test_indices,selected_features,with=F],feat_matrix[train_indices,selected_features,with=F])
	ranked = t(apply(dist_m,MARGIN=1,function(r)rank(r,ties.method="min")))
	rownames(ranked)=test_indices
	colnames(ranked)=train_indices
	rownames(dist_m)=test_indices
	colnames(dist_m)=train_indices
	ranked_dist=data.table(cbind(melt(dist_m),rank=as.vector(ranked)))
	setnames(ranked_dist,c("test_i","train_i","dist",'rank'))
	ranked_dist=ranked_dist[rank<=n_closer]
	if(!is.null(selected_labels)){

		test_attributes=cbind(feat_matrix[test_indices,selected_labels,with=F],i=test_indices);
		train_attributes=cbind(feat_matrix[train_indices,selected_labels,with=F],i=train_indices);
		setnames(test_attributes,paste("test",colnames(test_attributes),sep="_"));
		setnames(train_attributes,paste("train",colnames(train_attributes),sep="_"));

		ranked_dist=merge(merge(ranked_dist,test_attributes,by='test_i'),train_attributes,by='train_i')
	}
	return(ranked_dist)

}

{ # data loading 
input_pattern="1111"
mixed_3mers = fread(paste('../140613_more_patterns/spaced_kmers/long_kmers_mVar_n30.fa_',input_pattern,".csv",sep=""))

mixed_3mers = dcast.data.table(mixed_3mers,pattern+GC+species+sequence_description~kmer,value.var='count',fill=0)

seq_attributes=limma::strsplit2(mixed_3mers$sequence_description,"_")
mixed_3mers$class=factor(seq_attributes[,2])
mixed_3mers$species=factor(seq_attributes[,1])
mixed_3mers$pattern=factor(mixed_3mers$pattern)

mixed_3mers$GC_f=cut(mixed_3mers$GC, quantile(mixed_3mers$GC,prob=seq(0,1,0.1)),include.lowest=T)


setkey(mixed_3mers,"sequence_description")

mixed_3mers_label_cols=c('GC_f','class',"pattern",'GC','species','sequence_description')
mixed_3mers_count_cols=colnames(mixed_3mers)
mixed_3mers_count_cols=setdiff(mixed_3mers_count_cols,mixed_3mers_label_cols)
}


test_res=compute_closest_train_test(mixed_3mers,train=1:nrow(mixed_3mers),test=1:5,mixed_3mers_count_cols,5,mixed_3mers_label_cols)
test_res=compute_closest_train_test(mixed_3mers,train=1:nrow(mixed_3mers)/2,test=nrow(mixed_3mers)/2:nrow(mixed_3mers),mixed_3mers_count_cols,5,mixed_3mers_label_cols)


## We simulate a CV setup 


all_contigs=unique(mixed_3mers[,list(sequence_description,class)])

all_res=data.table()

training_ratio=0.9

splits=createDataPartition(all_contigs$class,p=training_ratio,times=5)
# for each split, compute the sparse distance matrix 


for( split_index in 1:length(splits)){
	cat(split_index,"\n")
	split1_training_idx=splits[[split_index]]
	split1_testing_idx=setdiff(1:nrow(all_contigs),split1_training_idx)
	n_closer_neighbors=60

	input_size=length(split1_testing_idx)
	nas_to_add=((input_size %/% n_cores)+1)*n_cores - input_size
	input_splits= matrix(c(split1_testing_idx,rep(NA,nas_to_add)),byrow=T,ncol=n_cores)

	res= data.table(foreach(i=iter(input_splits,by='col'),.combine='rbind',.inorder=F) %dopar% compute_closest_train_test(mixed_3mers,test=i[,1],train=split1_training_idx,mixed_3mers_count_cols,n_closer=n_closer_neighbors,selected_labels=c("class","species","sequence_description")))

	res$metric=factor('euclidean',levels=c('euclidean','length_norm','scaled'))
	res$pattern=input_pattern
	res$training_ratio=training_ratio
	res$split_index=split_index
	all_res=rbind(all_res,res)
}

# we classify using 11 neighbors 
res_classified=all_res[rank<=11,list(predicted_class=names(sort(table(train_class))[4])),by=list(test_sequence_description,test_class,split_index,training_ratio)]
by_splits_train=by(res_classified,list(res_classified$split_index,res_classified$training_ratio),function(dt){as.data.frame(t(confusionMatrix(dt$predicted_class,dt$test_class)$overall))})
by_splits_train=cbind(expand.grid('idx'=rownames(by_splits_train),'ratio'=colnames(by_splits_train)),do.call('rbind',by_splits_train))


# We plot as a function of the number of neihbors 


by_n_neigh=data.table(ldply(1:50,function(a_rank){
	res_classified=all_res[rank<=a_rank,list(predicted_class=names(sort(table(train_class))[4])),by=list(test_sequence_description,test_class,split_index,training_ratio)]
	by_splits_train=by(res_classified,list(res_classified$split_index,res_classified$training_ratio),function(dt){as.data.frame(t(confusionMatrix(dt$predicted_class,dt$test_class)$overall))})
	by_splits_train=cbind(expand.grid('idx'=rownames(by_splits_train),'ratio'=colnames(by_splits_train)),do.call('rbind',by_splits_train))

	by_splits_train$n_neighbors=a_rank
	return(by_splits_train)
},.parallel=T))
ggplot(by_n_neigh,aes(x=n_neighbors,y=Kappa,colour=idx))+geom_point()+geom_line()

## What's the timing for a single core single instance ? How many subset of features can we scan ? 
## What if we split by species ? 
## What if we split by contig ? 
# What if we remove some kMer cols ? 
# Any contig causing a lot of missclassification? (short one)
# Normalization helps ?
# Could spectral clust help ? 
# What about other distances from tomovic ? 
