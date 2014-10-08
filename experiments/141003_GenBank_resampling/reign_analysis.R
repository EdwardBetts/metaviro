### hack hack hack hack hack hack hack 
### Since we the target classification is used in formulas, we cannot pass it as an argument :( 
### However we bypass this problem by creating a proxy variable with a fixed name; having the values of the classification variable of interest ! 


normalized_kmers$classLabels=normalized_kmers[,classLabelsTag,with=F]
all_permutations_classes=data.table()
all_permutations_classes_means=data.table()

all_classes_kDN=data.table()
all_classes_kDN_means=data.table()

labels_to_return=c("gi","domain")


most_frequent_classes=tail(normalized_kmers[grep("unclassified|unassigned",classLabels,invert=T)][,.N,by=classLabels][order(N)],n=4)$classLabels
min_size=normalized_kmers[classLabels %in% most_frequent_classes,.N,by=classLabels][,min(N)]

for(group_sampling_id in 1:N_GROUP_SAMPLING){
	loginfo("[%s]Doing downsampling %d/%d",classLabelsTag,group_sampling_id,N_GROUP_SAMPLING)

	downsample_idx=unlist(llply(most_frequent_classes,function(g){sample(normalized_kmers[classLabels == g,idx],min_size)}))
	normalized_kmers_ds_SG=normalized_kmers[downsample_idx]

	input_size=nrow(normalized_kmers_ds_SG)
	nas_to_add=((input_size %/% n_cores)+1)*n_cores - input_size
	input_splits= matrix(c(1:input_size,rep(NA,nas_to_add)),byrow=T,ncol=n_cores)

	res_maj_sub_group= data.table(foreach(i=iter(input_splits,by='col'),.combine='rbind',.inorder=F) %dopar% 
		compute_closest2(normalized_kmers_ds_SG,kmers_columns,i[,1],N_NEAREST_NEIGHBORS,selected_labels=c(labels_to_return,"classLabels")))
	loginfo("[%s]Doing downsampling %d/%d: Computed nearest neighbors for kDN",classLabelsTag,group_sampling_id,N_GROUP_SAMPLING)



	# Tabulate 
	res_maj_sub_group_no_diag=res_maj_sub_group[src!=tgt & rank<=N_NEAREST_NEIGHBORS]

	res_maj_sub_group_neighbors_group=data.table(dcast(res_maj_sub_group_no_diag,src+src_gi+src_classLabels~tgt_classLabels,fun.aggregate =length,value.var="tgt_gi"))
	res_maj_sub_group_neighbors_group$src_classLabels=as.character(res_maj_sub_group_neighbors_group$src_classLabels)

	res_maj_sub_group_neighbors_sub_group_m=melt(res_maj_sub_group_neighbors_group,id.vars=c("src","src_gi","src_classLabels"))
	res_maj_sub_group_neighbors_sub_group_m[,mismatch:=F]
	res_maj_sub_group_neighbors_sub_group_m[src_classLabels!=variable,mismatch:=T]

	res_maj_sub_group_neighbors_sub_group_matches=res_maj_sub_group_neighbors_sub_group_m[,sum(value),by=list(src,src_gi,src_classLabels,variable,mismatch)]
	res_maj_sub_group_neighbors_sub_group_matches_disagreeing=res_maj_sub_group_neighbors_sub_group_matches[mismatch==T]

	sub_groups_summaries=res_maj_sub_group_neighbors_sub_group_matches_disagreeing[,list("kDN"=sum(V1)),by=list(src,src_gi,src_classLabels,mismatch)]
	sub_groups_summaries$group_sampling_id=group_sampling_id
	all_classes_kDN=rbind(all_classes_kDN,sub_groups_summaries)
	
	sub_group_means=sub_groups_summaries %>% group_by(src_classLabels) %>% dplyr::summarize(mean_kDN=mean(kDN),counts=n(),sd_kDN=sd(kDN),median_kDN=median(kDN))
	sub_group_means$group_sampling_id=group_sampling_id
	all_classes_kDN_means=rbind(all_classes_kDN_means,sub_group_means)







	# do the permutation analysis 

	for (perm_id in 1:N_PERMUTATIONS){
		loginfo("[%s]Doing permutation %d/%d",classLabelsTag,perm_id,N_PERMUTATIONS)

		res_maj_by_sub_group_perm=res_maj_sub_group[src!=tgt & rank<=N_NEAREST_NEIGHBORS]
		res_maj_by_sub_group_perm$tgt_classLabels=sample(res_maj_sub_group[src!=tgt & rank<=N_NEAREST_NEIGHBORS,tgt_classLabels])

		res_by_sub_group_perm=data.table(dcast(res_maj_by_sub_group_perm,src+src_gi+src_classLabels~tgt_classLabels,fun.aggregate =length,value.var="tgt_gi"))
		res_by_sub_group_perm$src_classLabels=as.character(res_by_sub_group_perm$src_classLabels)
		res_by_sub_group_perm[,.N,by=src_classLabels]


		res_by_sub_group_perm_m=melt(res_by_sub_group_perm,id.vars=c("src","src_gi","src_classLabels"))
		res_by_sub_group_perm_m[,mismatch:=F]
		res_by_sub_group_perm_m[src_classLabels!=variable,mismatch:=T]

		res_by_sub_group_perm_matches=res_by_sub_group_perm_m[,sum(value),by=list(src,src_gi,src_classLabels,variable,mismatch)]
		res_by_sub_group_perm_matches_disagreeing=res_by_sub_group_perm_matches[mismatch==T]

		classLabelss_summaries_perm=res_by_sub_group_perm_matches_disagreeing[,list("kDN"=sum(V1)),by=list(src,src_gi,src_classLabels,mismatch)]
		classLabelss_summaries_perm$perm_id=perm_id


		sub_group_means_perm=classLabelss_summaries_perm %>% group_by(src_classLabels) %>% dplyr::summarize(mean_kDN=mean(kDN),counts=n(),sd_kDN=sd(kDN),median_kDN=median(kDN))
		sub_group_means_perm$perm_id=perm_id
		
		all_permutations_classes=rbind(all_permutations_classes,classLabelss_summaries_perm)
		all_permutations_classes_means=rbind(all_permutations_classes_means,sub_group_means_perm)
	}


}

# combined_means=rbind(data.table(all_permutations_classes_means,model="null",group_sampling_id=-1), data.table(all_classes_kDN_means,model="observed",perm_id=-1),use.names=T)

# pdf(paste(DOMAIN,classLabelsTag,"nulls_kDN.pdf",sep="_"),w=16,h=10)
# g=ggplot(combined_means,aes(x=src_classLabels,y=mean_kDN,colour=model,size=model))+geom_point()+stat_summary(data=combined_means[model=="null"],fun.data = "mean_cl_boot", colour = "red",size=1)+ggtitle("Distribution of mean observed kDN values vs mean kDN values sampled under the null")+scale_x_discrete(name=classLabelsTag)+scale_y_continuous(name="mean kDN",limits=c(0,N_NEAREST_NEIGHBORS))+scale_colour_manual(values=c("grey40","black"))+scale_size_manual(values=c(1,3))
# print(g)
# foo=dev.off()

# pdf(paste(DOMAIN,classLabelsTag,"distribution_kDN.pdf",sep="_"),w=16,h=10)
# g=ggplot(all_classes_kDN,aes(x=factor(kDN),fill=src_classLabels))+geom_bar()+scale_x_discrete(limits=factor(0:N_NEAREST_NEIGHBORS),breaks=factor(seq(0,N_NEAREST_NEIGHBORS,10)),name="kDN")+scale_fill_discrete(name=classLabelsTag)+ggtitle(paste("Distribution of kDN values over all",N_GROUP_SAMPLING," downsampling"))
# print(g)
# foo=dev.off()


# # median_kDN_over_resampling=all_classes_kDN[,list("kDN"=as.numeric(median(kDN))),by=list(src_gi,src_classLabels)]

# # ggplot(median_kDN_over_resampling,aes(x=factor(kDN),fill=src_classLabels))+geom_bar()+scale_x_discrete(limits=factor(0:N_NEAREST_NEIGHBORS),breaks=factor(seq(0,N_NEAREST_NEIGHBORS,10)),name="kDN")+scale_fill_discrete(name=classLabelsTag)+ggtitle("Distribution of median kDN values per species over all downsampling")
# # print(g)
# # foo=dev.off()

# # ggplot(all_classes_kDN[src_gi %in% all_classes_kDN[,.N,by=src_gi][N>=10,src_gi]],aes(x=src_gi,y=kDN))+geom_point()+facet_wrap(~src_classLabels,scale="free_x")

# # # Second output 
# # pdf(paste0(DOMAIN,"_sub_group_kDN.pdf"),w=24,h=16)
# # # ggplot(groups_summaries,aes(x=factor(kDN,levels=0:N_NEAREST_NEIGHBORS),fill=))+geom_histogram()
# # g=
# # print(g)
# # foo=dev.off()



# # # Save the kDN results 

# write.csv(combined_means,file=paste(DOMAIN,classLabelsTag,"mean_kDN.csv",sep="_"),row.names=F)
# write.csv(all_classes_kDN,file=paste(DOMAIN,classLabelsTag,"by_gi_kDN.csv",sep="_"),row.names=F)
