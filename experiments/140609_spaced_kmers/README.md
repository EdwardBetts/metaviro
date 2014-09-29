SPACED (cf [spaced]) generates patterns for k-mers masking based on non-periodic 

We hypothesize that such spaced k-mers generate distributions matching tax assignment 

* optimal length 
The optimal length of k-mers are determined empirically using mitochondrial genomes, and was found to be k==9 or k==8 (synthetic examples)

* Number of pattern 

User-defined 

For 16k nt sequences, between 60 and 70 patterns showed best results (as measured by phylogenetic trees reconstructions)

* Pattern generation 

Patterns are of length between 9 and 39 with as much 9 1's 
Distribution is uniform ? Might be 


As many 1s as in the optimal pattern (k==9) thus as much 
Actually we always the same number of 1's
But the number of 0's is random, such that it's less or equal than 1's 

Position of 1's and 0's are random 


* Length distribution 
K-mers of length up to 39 were used in their studies 
For length > 20, eucl distance yield best results 

* Pattern weights
No idea 

* Normalization of pattern counts 
No idea 

* Multiple vs single patterns 

Multiple is better in almost all cases
(Single <=> only one pattern)

* Vs classical approaches 

Usually, people do either pair-wise alignment or multiple alignment between the input sequences. 
Observed that MLikelihood+ multiple alignment performed better than multiple pattern + eucl distance 

* Protein alignment vs nucleotide alignment 
For protein alignment k==4  was used, with a max length of 34 


* [spaced] 