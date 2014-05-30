METAVIROHOME:=~/metaviro
BIN:=~/metaviro/bin
DATA:=$(METAVIROHOME)/data

full_ncbi_contigs_4k: $(DATA)/full_ncbi/archea.fa $(DATA)/full_ncbi/bact.fa $(DATA)/full_ncbi/euk.fa $(DATA)/full_ncbi/viruses.fa 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/archea.fa -ARs -k "archea" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/bact.fa -ARs -k "bact" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/euk.fa -ARs -k "euk" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/viruses.fa -ARs -k "viruses" -n 4000 -o $@ 

full_ncbi_contigs_4k_3mers.csv: full_ncbi_contigs_4k
	$(BIN)/kmerize.py -k3 $< -o $@