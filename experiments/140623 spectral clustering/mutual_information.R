# we want to determine if there are kMers with stronger mutual information with the class labels than the other 


{ #setup
	library(caret)
	library(gridExtra)
	library(fields)
	library(party)
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



## normalize 

normalized_counts=mixed_3mers[,mixed_3mers_count_cols,with=F] / rowSums(mixed_3mers[,mixed_3mers_count_cols,with=F])
normalized_counts$class=mixed_3mers$class

## Correlation 

kmer_cor_with_class=data.table(t(cor(as.numeric(normalized_counts[,class]),as.matrix(normalized_counts[,mixed_3mers_count_cols,with=F]))),keep.rownames=T)

kmer_cor_with_class[order(V1)]
top10=kmer_cor_with_class[order(V1),rn] %.% tail(n=10)



ggplot(melt(normalized_counts[,c(top10,"class"),with=F]),aes(x=class,y=value,colour=class))+geom_point()+facet_wrap(~variable)

## rforest on the top 5 
top5=kmer_cor_with_class[order(V1),rn] %.% tail(n=5)
t=ctree(data=normalized_counts[,c(top5,"class"),with=F],class~.)

## We make an eikosogram 
