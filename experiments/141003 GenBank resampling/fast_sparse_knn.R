library(data.table,quietly=T,warn.conflicts = F)
library(fields,quietly=T,warn.conflicts = F)
library(reshape2,quietly=T,warn.conflicts = F)

# we assume that the feat matrix is a data table and has an id column 
compute_closest2 = function(feat_matrix,selected_features,selected_obs,n_closer,selected_labels=NA){

	src_idx=na.omit(selected_obs)
	dist_m=rdist(feat_matrix[src_idx,selected_features,with=F],feat_matrix[,selected_features,with=F])
	rownames(dist_m)=src_idx
	ranked=data.table(melt(apply(dist_m,MARGIN=1,function(r)rank(r,ties.method="min"))))
	setnames(ranked,c("tgt","src","rank"))
	unique(ranked$src)
	setdiff(ranked$src,src_idx)
	setdiff(src_idx,ranked$src)

	ranked_filt=ranked[rank<=n_closer]
	setdiff(ranked_filt$src,src_idx)
	setdiff(src_idx,ranked_filt$src)

	unique(feat_matrix[,c("idx",selected_labels),with=F])

	src_attributes=feat_matrix[ranked_filt$src,c("idx",selected_labels),with=F]
	tgt_attributes=feat_matrix[ranked_filt$tgt,c("idx",selected_labels),with=F]

	setnames(src_attributes,selected_labels,paste("src",selected_labels,sep="_"));
	setnames(tgt_attributes,selected_labels,paste("tgt",selected_labels,sep="_"));
	setnames(src_attributes,"idx","src")
	setnames(tgt_attributes,"idx","tgt")


	ranked_filt=cbind(ranked_filt,src_attributes,tgt_attributes)
	return(ranked_filt)
}



determine_disagreeing = function(n_closer_neighbors){
	most_frequent_groups=tail(normalized_kmers[,.N,by=Group][order(N)],n=4)$Group
	min_size=min(tail(normalized_kmers[,.N,by=Group][order(N)],n=4)$N)
	labels_to_return=c("Group","SubGroup","Host","Status","Segments","Genes","Proteins","org","org.GC","species")


	downsample_idx=unlist(llply(most_frequent_groups,function(g){sample(normalized_kmers[Group == g,idx],min_size)}))
	normalized_kmers_ds_SG=normalized_kmers[downsample_idx]
	normalized_kmers_ds_SG[,.N,by=Group]
	normalized_kmers_ds_SG$idx=1:nrow(normalized_kmers_ds_SG)
	input_size=nrow(normalized_kmers_ds_SG)
	nas_to_add=((input_size %/% n_cores)+1)*n_cores - input_size
	input_splits= matrix(c(1:input_size,rep(NA,nas_to_add)),byrow=T,ncol=n_cores)
	res_maj_by_group= data.table(foreach(i=iter(input_splits,by='col'),.combine='rbind',.inorder=F) %dopar% compute_closest2(normalized_kmers_ds_SG,kmers_columns,i[,1],n_closer_neighbors,selected_labels=labels_to_return))


	# permute class labels 
	res_maj_by_group_orig=res_maj_by_group
	res_maj_by_group$tgt_Group=sample(res_maj_by_group$tgt_Group)
	# unpermute 
	# res_maj_by_group$tgt_Group=res_maj_by_group_orig$tgt_Group


	res_maj_by_group_no_diag=res_maj_by_group[src!=tgt & rank<=n_closer_neighbors]
	# res_maj_by_group_no_diag[src==1]

	res_maj_by_group_neighbors_group=data.table(dcast(res_maj_by_group_no_diag,src+src_Group~tgt_Group))
	res_maj_by_group_neighbors_group$src_Group=as.character(res_maj_by_group_neighbors_group$src_Group)
	res_maj_by_group_neighbors_group[,.N,by=src_Group]


	res_maj_by_group_neighbors_sub_group_m=melt(res_maj_by_group_neighbors_group,id.vars=c("src","src_Group"))
	res_maj_by_group_neighbors_sub_group_m[,mismatch:=F]
	res_maj_by_group_neighbors_sub_group_m[src_Group!=variable,mismatch:=T]

	res_maj_by_group_neighbors_sub_group_matches=res_maj_by_group_neighbors_sub_group_m[,sum(value),by=list(src,src_Group,variable,mismatch)]
	res_maj_by_group_neighbors_sub_group_matches_disagreeing=res_maj_by_group_neighbors_sub_group_matches[mismatch==T]
	sub_groups_summaries=res_maj_by_group_neighbors_sub_group_matches_disagreeing[,sum(V1),by=list(src,src_Group,mismatch)]
	print(tbl_df(sub_groups_summaries) %.% group_by(src_Group) %.% summarise(a_mean=mean(V1)))
	return(res_maj_by_group_neighbors_sub_group_matches)
}
