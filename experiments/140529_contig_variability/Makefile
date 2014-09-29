full_ncbi_contigs_per_sequences.fa:
	../../bin/contig_sampling.py -n 10 -S ../../data/full_ncbi/archea.fa -ARs -k "arch" -o $@
	../../bin/contig_sampling.py -n 100 -S ../../data/full_ncbi/euk.fa -ARs -k "euk" -o $@ 
	../../bin/contig_sampling.py -n 10 -S ../../data/full_ncbi/bact.fa -ARs -k "bact" -o $@
	../../bin/contig_sampling.py -n 5 -S ../../data/full_ncbi/viruses.fa -ARs -k "virus" -o $@

full_ncbi_contigs_per_sequences_3mer.csv: full_ncbi_contigs_per_sequences.fa
	../../bin/kmerize.py -k3 $< -o $@ 
	
