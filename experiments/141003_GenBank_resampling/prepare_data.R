suppressMessages(library(reshape2,quietly = T,warn.conflicts = F))
suppressMessages(library(plyr,quietly = T,warn.conflicts = F))
suppressMessages(library(dplyr,quietly = T,warn.conflicts = F))
suppressMessages(library(ggplot2,quietly = T,warn.conflicts = F))
suppressMessages(library(data.table,quietly = T,warn.conflicts = F))
suppressMessages(library(foreach,quietly = T,warn.conflicts = F))
suppressMessages(library(doParallel,quietly = T,warn.conflicts = F))





registerDoParallel(cores=n_cores)
# We load the source data 

all_instances_kmers=fread(SOURCE_KMERS)

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
domain_kmers=all_instances_kmers[domain==DOMAIN] # 10814 contigs 
setkey(domain_kmers,"gi")

# # We load the GI to TaxID conversion 
load(GI_MAPPING)

# setkey(virus_annotations_gi,"gi")
annotated_kmers=merge(domain_kmers,GenBank_taxID2gi,by="gi")
annotated_kmers=merge(annotated_kmers,GenBank_annotations,by="TaxID")

# normalize and type 
normalized_kmers_counts=annotated_kmers[,kmers_columns,with=F]/annotated_kmers$sequence_length

normalized_kmers=cbind(normalized_kmers_counts, annotated_kmers[,list(gi,BioProject.ID,Group,SubGroup,Status,Genes,Proteins,"org"=Organism.Name,"org.GC"=GC.,"size"=Size)] )

normalized_kmers$org.GC=as.numeric(normalized_kmers$org.GC)
normalized_kmers$size=as.numeric(normalized_kmers$size)
normalized_kmers$Genes=as.numeric(normalized_kmers$Genes)
normalized_kmers$Proteins=as.numeric(normalized_kmers$Proteins)

normalized_kmers$Status=as.factor(normalized_kmers$Status)
normalized_kmers$Group=as.factor(normalized_kmers$Group)
normalized_kmers$SubGroup=as.factor(normalized_kmers$SubGroup)
normalized_kmers$gi=as.factor(normalized_kmers$gi)
normalized_kmers$org=as.factor(normalized_kmers$org)
normalized_kmers$idx=1:nrow(normalized_kmers)
