# On Centos, require 
# module load python/2.7.3
METAVIROHOME:=~/metaviro
BIN:=~/metaviro/bin
DATA:=$(METAVIROHOME)/data

.PHONY: all clean 
all: long_kmers.fa.csv
clean:
	rm long_kmers.fa long_kmers.fa.csv

long_kmers.fa: $(DATA)/full_ncbi/archea.fa $(DATA)/full_ncbi/bact.fa $(DATA)/full_ncbi/euk.fa $(DATA)/full_ncbi/viruses.fa 
	$(BIN)/contig_sampling.py -M 30 -S -n 20 -s -k "archea" ../../data/full_ncbi/archea.fa -o $@
	$(BIN)/contig_sampling.py -A -M 30 -S -n 20 -s -k "bact" ../../data/full_ncbi/bact.fa -o $@
	$(BIN)/contig_sampling.py -A -M 30 -S -n 20 -s -k "euk" ../../data/full_ncbi/euk.fa -o $@
	$(BIN)/contig_sampling.py -A -M 30 -S -n 20 -s -k "viruses" ../../data/full_ncbi/viruses.fa -o $@


long_kmers.fa.csv: long_kmers.fa
	$(BIN)/kmerize.py -k10 -L $< -o $@
	$(BIN)/kmerize.py -A -k11 -L $< -o $@
	$(BIN)/kmerize.py -A -k10 -L -S $< -o $@
	$(BIN)/kmerize.py -A -k11 -L -S $< -o $@	

long_kmers_b2.fa: $(DATA)/full_ncbi/archea.fa $(DATA)/full_ncbi/bact.fa $(DATA)/full_ncbi/euk.fa $(DATA)/full_ncbi/viruses.fa 
	$(BIN)/contig_sampling.py -M 50 -S -n 5 -s -k "archea" ../../data/full_ncbi/archea.fa -o $@
	$(BIN)/contig_sampling.py -A -M 50 -S -n 5 -s -k "bact" ../../data/full_ncbi/bact.fa -o $@
	$(BIN)/contig_sampling.py -A -M 50 -S -n 5 -s -k "euk" ../../data/full_ncbi/euk.fa -o $@
	$(BIN)/contig_sampling.py -A -M 50 -S -n 5 -s -k "viruses" ../../data/full_ncbi/viruses.fa -o $@


long_kmers_b2.fa.csv: long_kmers_b2.fa
	$(BIN)/kmerize.py -k10 -L $< -o $@
	$(BIN)/kmerize.py -A -k11 -L $< -o $@
	$(BIN)/kmerize.py -A -k10 -L -S $< -o $@
	$(BIN)/kmerize.py -A -k11 -L -S $< -o $@	

long_kmers_b3.fa: $(DATA)/full_ncbi/archea.fa $(DATA)/full_ncbi/bact.fa $(DATA)/full_ncbi/euk.fa $(DATA)/full_ncbi/viruses.fa 
	$(BIN)/contig_sampling.py -M 100 -S -n 4 -s -k "archea" ../../data/full_ncbi/archea.fa -o $@
	$(BIN)/contig_sampling.py -A -M 100 -S -n 4 -s -k "bact" ../../data/full_ncbi/bact.fa -o $@
	$(BIN)/contig_sampling.py -A -M 100 -S -n 4 -s -k "euk" ../../data/full_ncbi/euk.fa -o $@
	$(BIN)/contig_sampling.py -A -M 100 -S -n 4 -s -k "viruses" ../../data/full_ncbi/viruses.fa -o $@


long_kmers_b3.fa.csv: long_kmers_b3.fa
	$(BIN)/kmerize.py -k10 -L $< -o $@
	$(BIN)/kmerize.py -A -k11 -L $< -o $@
	$(BIN)/kmerize.py -A -k10 -L -S $< -o $@
	$(BIN)/kmerize.py -A -k11 -L -S $< -o $@	