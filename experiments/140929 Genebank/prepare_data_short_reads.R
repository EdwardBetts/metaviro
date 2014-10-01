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

short_virus_kmers=fread("../../data/Genebank/GB_virus_contigs_n_4_short_reads.csv")
# We build a matrix of normalized kmers 

kmers_columns=colnames(short_virus_kmers)[5:ncol(short_virus_kmers)]
sequence_attributes=limma::strsplit2(short_virus_kmers$sequence_description,"_",fixed=T)

short_virus_kmers$species=factor(sequence_attributes[,1])
short_virus_kmers$class=factor(sequence_attributes[,2])
short_virus_kmers$contig_index=as.numeric(sequence_attributes[,4])

viruses_bp_accession_numbers=cbind("species"=sequence_attributes[,1],data.table(limma::strsplit2(sequence_attributes[,5],"/|-")))
setnames(viruses_bp_accession_numbers,"V2","BioProject.ID")
viruses_bp_accession_numbers$BioProject.ID=as.integer(as.character(viruses_bp_accession_numbers$BioProject.ID))
viruses_bp_accession_numbers=unique(viruses_bp_accession_numbers[V1=="viruses"])
short_virus_kmers=merge(short_virus_kmers,viruses_bp_accession_numbers[,list(species,BioProject.ID)],by="species")
virus_annotations=fread("../../data/GENOME_REPORTS/viruses.txt")
setnames(virus_annotations,make.names(colnames(virus_annotations)))
# virus_annotations$BioProject.ID=as.integer(as.character(virus_annotations$BioProject.ID))

short_virus_kmers=merge(short_virus_kmers,virus_annotations,by="BioProject.ID")

# We remove classes with not enough instances 
short_virus_kmers[Host=="algea",Host:="algae"]
short_virus_kmers[Host=="human",Host:="vertebrates, human"]
groups_to_remove=short_virus_kmers[,.N,by=list(Group)][order(N)][N<=65,Group]
hosts_to_remove=short_virus_kmers[,.N,by="Host"][order(N)][N<=65,Host]

short_virus_kmers[Host %in% hosts_to_remove,Host:="NA"]
short_virus_kmers[Group %in% groups_to_remove,Group:="unclassified"]

# normalize and type 
short_virus_kmers_norm_counts=short_virus_kmers[,kmers_columns,with=F]/short_virus_kmers$sequence_length

short_virus_kmers_norm=cbind(short_virus_kmers_norm_counts, short_virus_kmers[,list(species,BioProject.ID,Group,SubGroup,Host,Status,"Segments"=Segmemts,Genes,Proteins,"org"=X.Organism.Name,"org.GC"=GC.,"size"=Size..Kb.)] )

short_virus_kmers_norm$org.GC=as.numeric(short_virus_kmers_norm$org.GC)
short_virus_kmers_norm$size=as.numeric(short_virus_kmers_norm$size)
short_virus_kmers_norm$Genes=as.numeric(short_virus_kmers_norm$Genes)
short_virus_kmers_norm$Proteins=as.numeric(short_virus_kmers_norm$Proteins)

short_virus_kmers_norm$Status=as.factor(short_virus_kmers_norm$Status)
short_virus_kmers_norm$Group=as.factor(short_virus_kmers_norm$Group)
short_virus_kmers_norm$SubGroup=as.factor(short_virus_kmers_norm$SubGroup)
short_virus_kmers_norm$Host=as.factor(short_virus_kmers_norm$Host)
short_virus_kmers_norm$species=as.factor(short_virus_kmers_norm$species)
short_virus_kmers_norm$org=as.factor(short_virus_kmers_norm$org)
short_virus_kmers_norm$idx=1:nrow(short_virus_kmers_norm)


# PCA 

pc_kmers=princomp(short_virus_kmers_norm_counts)
coord_pca=predict(pc_kmers)[,1:3]
short_virus_kmers_norm=cbind(short_virus_kmers_norm,coord_pca)
