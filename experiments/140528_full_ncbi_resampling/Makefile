# On Centos, require 
# module load python/2.7.3
METAVIROHOME:=~/metaviro
BIN:=~/metaviro/bin
DATA:=$(METAVIROHOME)/data

full_ncbi_contigs_4k.fa: $(DATA)/full_ncbi/archea.fa $(DATA)/full_ncbi/bact.fa $(DATA)/full_ncbi/euk.fa $(DATA)/full_ncbi/viruses.fa 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/archea.fa -ARsP -k "archea" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/bact.fa -ARsP -k "bact" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/euk.fa -ARsP -k "euk" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/viruses.fa -ARsP -k "viruses" -n 4000 -o $@ 


full_ncbi_contigs_4k_v2.fa: $(DATA)/full_ncbi/archea.fa $(DATA)/full_ncbi/bact.fa $(DATA)/full_ncbi/euk.fa $(DATA)/full_ncbi/viruses.fa 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/archea.fa -ARsP -k "archea" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/bact.fa -ARsP -k "bact" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/euk.fa -ARsP -k "euk" -n 4000 -o $@ 
	$(BIN)/contig_sampling.py $(DATA)/full_ncbi/viruses.fa -ARsP -k "viruses" -n 4000 -o $@ 


full_ncbi_contigs_4k_stats.txt: full_ncbi_contigs_4k.fa
	../../bin/fasta_statistics.py -i $< -d ";" -o $@
	R --slave --vanilla -e "library(data.table);stats=fread(\"$@\");summary(stats)"

full_ncbi_contigs_4k_v2_stats.txt: full_ncbi_contigs_4k_v2.fa
	../../bin/fasta_statistics.py -i $< -d ";" -o $@
	R --slave --vanilla -e "library(data.table);stats=fread(\"$@\");summary(stats)"

full_ncbi_contigs_4k_3mers.csv: full_ncbi_contigs_4k.fa
	$(BIN)/kmerize.py -k3 $< -o $@

full_ncbi_contigs_4k_v2_3mers.csv: full_ncbi_contigs_4k_v2.fa
	$(BIN)/kmerize.py -k3 $< -o $@


full_ncbi_contigs_4k_v2_7mers.csv: full_ncbi_contigs_4k_v2.fa
	$(BIN)/kmerize.py -k7 $< -o $@

full_ncbi_contigs_4k_v2_5mers.csv: full_ncbi_contigs_4k_v2.fa
	$(BIN)/kmerize.py -k5 $< -o $@