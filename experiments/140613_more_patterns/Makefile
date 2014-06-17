# On Centos, require 
# module load python/2.7.3
METAVIROHOME:=~/metaviro
BIN:=~/metaviro/bin
DATA:=$(METAVIROHOME)/data

.PHONY: all clean 
all: long_kmers.fa.csv
clean:
	rm long_kmers_m300_n4.fa long_kmers_m300_n4.fa.csv

long_kmers_m300_n4.fa: $(DATA)/full_ncbi/archea.fa $(DATA)/full_ncbi/bact.fa $(DATA)/full_ncbi/euk.fa $(DATA)/full_ncbi/viruses.fa 
	$(BIN)/contig_sampling.py -M 300 -S -n 4 -s -k "archea" ../../data/full_ncbi/archea.fa -o $@
	$(BIN)/contig_sampling.py -A -M 300 -S -n 4 -s -k "bact" ../../data/full_ncbi/bact.fa -o $@
	$(BIN)/contig_sampling.py -A -M 300 -S -n 4 -s -k "euk" ../../data/full_ncbi/euk.fa -o $@
	$(BIN)/contig_sampling.py -A -M 300 -S -n 4 -s -k "viruses" ../../data/full_ncbi/viruses.fa -o $@


long_kmers_m300_n4.fa.csv: long_kmers_m300_n4.fa
	$(BIN)/kmerize.py -k10 -L $< -o $@
	$(BIN)/kmerize.py -A -k11 -L $< -o $@
	$(BIN)/kmerize.py -A -k10 -L -S $< -o $@
	$(BIN)/kmerize.py -A -k11 -L -S $< -o $@	