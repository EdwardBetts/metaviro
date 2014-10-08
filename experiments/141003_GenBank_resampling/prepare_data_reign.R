
suppressMessages(library(reshape2,quietly = T,warn.conflicts = F))
suppressMessages(library(plyr,quietly = T,warn.conflicts = F))
suppressMessages(library(dplyr,quietly = T,warn.conflicts = F))
suppressMessages(library(ggplot2,quietly = T,warn.conflicts = F))
suppressMessages(library(data.table,quietly = T,warn.conflicts = F))
suppressMessages(library(foreach,quietly = T,warn.conflicts = F))
suppressMessages(library(doParallel,quietly = T,warn.conflicts = F))



n_cores=1
N_PERMUTATIONS=10
N_GROUP_SAMPLING=1
N_NEAREST_NEIGHBORS=73

N_MAJORITY_CLASS_LABELS=4

library(logging,quietly=T)
logReset()
addHandler(writeToConsole)
classLabelsTag="domain"
source("../fast_sparse_knn.R")
registerDoParallel(cores=n_cores)
# We load the source data 

all_instances_kmers=rbind(
	fread("archaea_k3_mers.tsv"),
	fread("bact_k3_mers.tsv"),
	fread("euk_k3_mers.tsv"),
	fread("viruses_k3_mers.tsv"))

# We build a matrix of normalized kmers 

kmers_columns=colnames(all_instances_kmers)[5:ncol(all_instances_kmers)]

sequence_attributes=data.table(cbind(
	limma::strsplit2(all_instances_kmers$path,"_|/",fixed=F),
	limma::strsplit2(all_instances_kmers$sequence_description,"|",fixed=T)))

sequence_attributes[,contig_index:=1:.N,by=V4]


all_instances_kmers$gi=as.integer(as.character(sequence_attributes[,V7]))
all_instances_kmers$domain=factor(sequence_attributes[,V4])
all_instances_kmers$batch=factor(sequence_attributes[,V2])
all_instances_kmers$mean_length=factor(sequence_attributes[,V3])
all_instances_kmers$contig_index=sequence_attributes$contig_index

# We restrict to domain
domain_kmers=all_instances_kmers
setkey(domain_kmers,"gi")


normalized_kmers_counts=domain_kmers[,kmers_columns,with=F]/domain_kmers$sequence_length

normalized_kmers=cbind(normalized_kmers_counts, domain_kmers[,list(gi,domain)])
normalized_kmers$idx=1:nrow(normalized_kmers)

rm(all_instances_kmers)
rm(sequence_attributes)
rm(domain_kmers)