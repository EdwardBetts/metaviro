library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(data.table)
library(stringr)
library(corrplot)
library(caret)

library(foreach)
library(doParallel)
n_cores=4
registerDoParallel(cores=n_cores)

setwd("~/metaviro/experiments/140929 Genebank/")


# We load the source data 

all_instances_kmers=fread("../../data/Genebank/GB_20mb_per_domain_k3.csv")
# We build a matrix of normalized kmers 

kmers_columns=colnames(all_instances_kmers)[5:ncol(all_instances_kmers)]
sequence_attributes=limma::strsplit2(all_instances_kmers$sequence_description,"_",fixed=T)

all_instances_kmers$species=factor(sequence_attributes[,1])
all_instances_kmers$class=factor(sequence_attributes[,2])
all_instances_kmers$contig_index=as.numeric(sequence_attributes[,4])

# We restrict to viruses and load the annotations
viruses_kmers=all_instances_kmers[class=="virus"] # 10814 contigs 
viruses_bp_accession_numbers=cbind("species"=sequence_attributes[,1],data.table(limma::strsplit2(sequence_attributes[,5],"/|-")))
setnames(viruses_bp_accession_numbers,"V2","BioProject.ID")
viruses_bp_accession_numbers$BioProject.ID=as.integer(as.character(viruses_bp_accession_numbers$BioProject.ID))
viruses_bp_accession_numbers=unique(viruses_bp_accession_numbers[V1=="viruses"])
viruses_kmers=merge(viruses_kmers,viruses_bp_accession_numbers[,list(species,BioProject.ID)],by="species")
virus_annotations=fread("../../data/GENOME_REPORTS/viruses.txt")
setnames(virus_annotations,make.names(colnames(virus_annotations)))
# virus_annotations$BioProject.ID=as.integer(as.character(virus_annotations$BioProject.ID))

annotated_virus_kmers=merge(viruses_kmers,virus_annotations,by="BioProject.ID")

# Any additional duplicates?
annotated_virus_kmers[,.N,by=species][order(N)] # 5539 species, with 1 or 2 contigs each 

# We remove classes with not enough instances 
annotated_virus_kmers[Host=="algea",Host:="algae"]
annotated_virus_kmers[Host=="human",Host:="vertebrates, human"]
groups_to_remove=annotated_virus_kmers[,.N,by=list(Group)][order(N)][N<=65,Group]
hosts_to_remove=annotated_virus_kmers[,.N,by="Host"][order(N)][N<=65,Host]
annotated_virus_kmers[Host %in% hosts_to_remove,Host:=NA]
annotated_virus_kmers[Group %in% groups_to_remove,Group:="unclassified"]

# normalize and type 
normalized_kmers_counts=annotated_virus_kmers[,kmers_columns,with=F]/annotated_virus_kmers$sequence_length

normalized_kmers=cbind(normalized_kmers_counts, annotated_virus_kmers[,list(species,BioProject.ID,Group,SubGroup,Host,Status,"Segments"=Segmemts,Genes,Proteins,"org"=X.Organism.Name,"org.GC"=GC.,"size"=Size..Kb.)] )

normalized_kmers$org.GC=as.numeric(normalized_kmers$org.GC)
normalized_kmers$size=as.numeric(normalized_kmers$size)
normalized_kmers$Genes=as.numeric(normalized_kmers$Genes)
normalized_kmers$Proteins=as.numeric(normalized_kmers$Proteins)

normalized_kmers$Status=as.factor(normalized_kmers$Status)
normalized_kmers$Group=as.factor(normalized_kmers$Group)
normalized_kmers$SubGroup=as.factor(normalized_kmers$SubGroup)
normalized_kmers$Host=as.factor(normalized_kmers$Host)
normalized_kmers$species=as.factor(normalized_kmers$species)
normalized_kmers$org=as.factor(normalized_kmers$org)
normalized_kmers$idx=1:nrow(normalized_kmers)


# PCA 

pc_kmers=princomp(normalized_kmers_counts)
coord_pca=predict(pc_kmers)[,1:3]
normalized_kmers=cbind(normalized_kmers,coord_pca)
