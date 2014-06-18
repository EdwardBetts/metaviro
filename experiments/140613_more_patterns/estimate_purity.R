library(data.table)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(reshape2)

setwd("~/Documents/metaviro/experiments/140613_more_patterns/")
# input_kmers=fread("long_kmers_m300_n4.fa_10101.csv")
# input_file="long_kmers_m300_n4.fa_10111001010100010011.csv"
# input_file="long_kmers_m300_n4.fa_11101011.csv"
# input_file="long_kmers_m300_n4.fa_11101011.csv"
# input_file="long_kmers_m300_n4.fa_11101111011.csv"
# input_file="long_kmers_m300_n4.fa_101111001.csv"
# input_file="long_kmers_m300_n4.fa_1010111001.csv"
# input_file="long_kmers_m300_n4.fa_10110000111.csv"
# input_file="long_kmers_m300_n4.fa_100001111001.csv"
# input_file="long_kmers_m300_n4.fa_111111111.csv"
# input_file="long_kmers_m300_n4.fa_11001111.csv"
# input_file="long_kmers_m300_n4.fa_10101.csv"
# input_file="long_kmers_m300_n4.fa_10110111.csv"
# input_file="long_kmers_m300_n4.fa_101001111011.csv"
# input_file="long_kmers_m300_n4.fa_111011100111.csv"
# input_file="long_kmers_m300_n4.fa_1101111100001.csv"
# input_file="long_kmers_m300_n4.fa_101010100001111.csv"
input_file="long_kmers_m300_n4.fa_111101000010101.csv"



# Extract pattern to be used as plot title
pattern_name=strsplit(input_file,"_|\\.")[[1]][6]
#
input_kmers=fread(input_file)



seq_attributes=limma::strsplit2(input_kmers$sequence_description,"_")
input_kmers$class=factor(seq_attributes[,2])
input_kmers$species=factor(seq_attributes[,1])
input_kmers$pattern=factor(input_kmers$pattern)
input_kmers$kmer=factor(input_kmers$kmer)

n_classes=input_kmers[,.N,by=list(class,species)][,.N,by=class]
n_classes_l=as.list(n_classes$N)
names(n_classes_l)=as.character(n_classes$class)


kmer_by_classes=data.table(dcast(input_kmers[,sum(count),by=list(kmer,class)][order(V1)][V1>=20],kmer~class,value.var="V1"))
setkey(kmer_by_classes,"kmer")
setkey(input_kmers,"kmer")
long_kmer_retrieval=input_kmers[,length(unique(species)),by=list(kmer,class)][order(V1)]



long_kmer_retrieval_d = data.table(dcast(long_kmer_retrieval,kmer~class,value.var="V1",fill=0))

long_kmer_retrieval_d[order(viruses)]
long_kmer_retrieval_d[,virus_purity:=viruses/(archea+bact+euk+viruses)]
long_kmer_retrieval_d[,euk_purity:=euk/(archea+bact+euk+viruses)]
long_kmer_retrieval_d[,bact_purity:=bact/(archea+bact+euk+viruses)]
long_kmer_retrieval_d[,arch_purity:=archea/(archea+bact+euk+viruses)]

long_kmer_retrieval_d[,virus_recall:=viruses/n_classes_l$viruses]
long_kmer_retrieval_d[,virus_precision:=virus_purity]
long_kmer_retrieval_d[,virus_f1_score:=virus_recall*virus_precision/(virus_precision+virus_recall)]

long_kmer_retrieval_d[,euk_recall:=euk/n_classes_l$euk]
long_kmer_retrieval_d[,euk_precision:=euk_purity]
long_kmer_retrieval_d[,euk_f1_score:=euk_recall*euk_precision/(euk_precision+euk_recall)]

long_kmer_retrieval_d[,bact_recall:=bact/n_classes_l$bact]
long_kmer_retrieval_d[,bact_precision:=bact_purity]
long_kmer_retrieval_d[,bact_f1_score:=bact_recall*bact_precision/(bact_precision+bact_recall)]

long_kmer_retrieval_d[,arch_recall:=archea/n_classes_l$archea]
long_kmer_retrieval_d[,arch_precision:=arch_purity]
long_kmer_retrieval_d[,arch_f1_score:=arch_recall*arch_precision/(arch_precision+arch_recall)]



long_kmer_retrieval_d[viruses>=250][order(virus_purity)]
long_kmer_retrieval_d[!is.na(virus_f1_score)][order(virus_f1_score)]

write.table(long_kmer_retrieval_d,file=paste("processed/",input_file,"_retrieval_stats.csv",sep=""),quote=FALSE,row.names=F)

# Prepare a table with up to 4 champions (f1 score) per class
champions=unique(rbind(
	head(long_kmer_retrieval_d[!is.na(arch_f1_score)][order(arch_f1_score,decreasing=T)],n=2) %.% mutate("best"="A"),
	head(long_kmer_retrieval_d[!is.na(bact_f1_score)][order(bact_f1_score,decreasing=T)],n=2) %.% mutate("best"="B"),
	head(long_kmer_retrieval_d[!is.na(euk_f1_score)][order(euk_f1_score,decreasing=T)],n=2) %.% mutate("best"="E"),
	head(long_kmer_retrieval_d[!is.na(virus_f1_score)][order(virus_f1_score,decreasing=T)],n=2) %.% mutate("best"="V")
	))

annotations=geom_text(data=champions,aes(label=paste(kmer,best,sep="-")),vjust=1,hjust=-0.1,colour="firebrick4")

g1=ggplot(long_kmer_retrieval_d[viruses>=6],aes(x=viruses,y=virus_purity,size=virus_f1_score,color=euk_f1_score))+geom_point()+xlim(0,n_classes_l$viruses)+ylim(0,1)+scale_size_continuous(range=c(2,6))+annotations

g2=ggplot(long_kmer_retrieval_d[euk>=6],aes(x=euk,y=euk_purity,color=virus_f1_score,size=euk_f1_score))+geom_point()+xlim(0,n_classes_l$euk)+ylim(0,1)+scale_size_continuous(range=c(2,6))+annotations

g3=ggplot(long_kmer_retrieval_d[bact>=6],aes(x=bact,y=bact_purity,size=bact_f1_score,colour=virus_f1_score))+geom_point()+xlim(0,n_classes_l$bact)+ylim(0,1)+scale_size_continuous(range=c(2,6))+annotations

g4=ggplot(long_kmer_retrieval_d[archea>=6],aes(x=archea,y=arch_purity,size=arch_purity,colour=virus_f1_score))+geom_point()+xlim(0,n_classes_l$arc)+ylim(0,1)+scale_size_continuous(range=c(2,6))+annotations

pdf(file=paste("processed\\",input_file,"_retrieval_stats.pdf",sep=""),w=24,h=16)
grid.arrange(g1,g2,g3,g4,main=pattern_name)
dev.off()


## Prec recall plot 


g1=ggplot(long_kmer_retrieval_d[viruses>=6],aes(x=virus_precision,y=virus_recall,size=virus_f1_score,color=euk_f1_score))+geom_point()+xlim(0,1)+ylim(0,1)+scale_size_continuous(range=c(2,6))+annotations
g2=ggplot(long_kmer_retrieval_d[euk>=6],aes(x=euk_precision,y=euk_recall,size=euk_f1_score,color=virus_f1_score))+geom_point()+xlim(0,1)+ylim(0,1)+scale_size_continuous(range=c(2,6))+annotations
g3=ggplot(long_kmer_retrieval_d[bact>=6],aes(x=bact_precision,y=bact_recall,size=bact_f1_score,color=virus_f1_score))+geom_point()+xlim(0,1)+ylim(0,1)+scale_size_continuous(range=c(2,6))+annotations
g4=ggplot(long_kmer_retrieval_d[archea>=6],aes(x=arch_precision,y=arch_recall,size=arch_f1_score,color=euk_f1_score))+geom_point()+xlim(0,1)+ylim(0,1)+scale_size_continuous(range=c(2,6))+annotations

pdf(file=paste("processed\\",input_file,"_PR_retrieval_stats.pdf",sep=""),w=24,h=16)
grid.arrange(g1,g2,g3,g4,main=pattern_name)
dev.off()







