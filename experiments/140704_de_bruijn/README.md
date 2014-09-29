140704 De bruijn graph based classification 

We test a classifier based on graph features 
We want to determine whether a graph centric measure (associated with each k-mer) might be more appropriate than raw k-mer count .
When we build the De Bruin graph with 5-mers, we manage to reach a Kappa of 0.5 when we classify with 15-NN over the matrix of degree of the k-mers as features. 


	In [47]: confusion_matrix(testing_data["class"],predicted_classes)
	Out[47]: 
	array([[1330,  393,  672,   65],
	       [ 212, 2685,  505,   19],
	       [ 338,  460, 2808,   77],
	       [ 355,  704,  633,  741]])

	In [48]: class_map=dict(zip(set(testing_data["class"]),range(0,4)))

	In [49]: kappa([class_map[x] for x in testing_data["class"]],[class_map[x] for x in predicted_classes])
	Out[49]: 0.49212133306088668



	In [52]: scipy.stats.itemfreq(predicted_classes)
	Out[52]: 
	array([['archea', 2235],
	       ['bact', 4242],
	       ['euk', 4618],
	       ['viruses', 902]], dtype=object)

	In [53]: scipy.stats.itemfreq(testin)                                                                    
	testing_data  testing_set   

	In [53]: scipy.stats.itemfreq(testing_data["class"])
	Out[53]: 
	array([['archea', 2460],
	       ['bact', 3421],
	       ['euk', 3683],
	       ['viruses', 2433]], dtype=object)




When we use 7-mers, kappa is at 0 (everythin falls into the euk class)

The distribution of degree 