Analysis of PCA results on 3-mers showed a very strong GC effect captured by the first dimensions. 

* Can we correct for it? 
* Can we improve the classification performance (approx kappa==0.65 for 3-mers with 15 neighbors kNN) by accounting for it? 



# 15-NN on the 20 first PC vs full matrix of 4-mers


Performances with 0.8 training on the 20 PCA:

	> predicted_classes=predict(knn_performances,testExpr)
	> confusionMatrix(predicted_classes,testClass)
	Confusion Matrix and Statistics

	          Reference
	Prediction archea bact euk viruses
	   archea     403   43  72      47
	   bact        94  625  57     102
	   euk         93   48 678     107
	   viruses     38   38  82     358

	Overall Statistics
	                                          
	               Accuracy : 0.7154          
	                 95% CI : (0.6986, 0.7318)
	    No Information Rate : 0.3081          
	    P-Value [Acc > NIR] : < 2.2e-16       
	                                          
	                  Kappa : 0.6159          
	 Mcnemar's Test P-Value : 2.984e-10       

	Statistics by Class:

	                     Class: archea Class: bact Class: euk Class: viruses
	Sensitivity                 0.6417      0.8289     0.7627         0.5831
	Specificity                 0.9282      0.8813     0.8758         0.9304
	Pos Pred Value              0.7133      0.7118     0.7322         0.6938
	Neg Pred Value              0.9030      0.9357     0.8923         0.8919
	Prevalence                  0.2177      0.2614     0.3081         0.2128
	Detection Rate              0.1397      0.2166     0.2350         0.1241
	Detection Prevalence        0.1958      0.3043     0.3210         0.1789
	Balanced Accuracy           0.7850      0.8551     0.8192         0.7567

If we drop the first PC: 

	> confusionMatrix(predicted_classes,testClass)
	Confusion Matrix and Statistics

	          Reference
	Prediction archea bact euk viruses
	   archea     384   50  73      52
	   bact       101  610  82     105
	   euk         95   64 659     132
	   viruses     48   30  75     325

	Overall Statistics
	                                          
	               Accuracy : 0.6856          
	                 95% CI : (0.6683, 0.7025)
	    No Information Rate : 0.3081          
	    P-Value [Acc > NIR] : < 2.2e-16       
	                                          
	                  Kappa : 0.575           
	 Mcnemar's Test P-Value : 3.841e-15       

	Statistics by Class:

	                     Class: archea Class: bact Class: euk Class: viruses
	Sensitivity                 0.6115      0.8090     0.7413         0.5293
	Specificity                 0.9225      0.8649     0.8542         0.9326
	Pos Pred Value              0.6869      0.6793     0.6937         0.6799
	Neg Pred Value              0.8951      0.9275     0.8811         0.8799
	Prevalence                  0.2177      0.2614     0.3081         0.2128
	Detection Rate              0.1331      0.2114     0.2284         0.1127
	Detection Prevalence        0.1938      0.3113     0.3293         0.1657
	Balanced Accuracy           0.7670      0.8369     0.7977         0.7310



On the first 64 dimensions: 


	> confusionMatrix(predicted_classes,testClass)
	Confusion Matrix and Statistics

	          Reference
	Prediction archea bact euk viruses
	   archea     402   32  53      47
	   bact        93  643  74     104
	   euk         99   46 689     114
	   viruses     34   33  73     349

	Overall Statistics
	                                          
	               Accuracy : 0.722           
	                 95% CI : (0.7053, 0.7383)
	    No Information Rate : 0.3081          
	    P-Value [Acc > NIR] : < 2.2e-16       
	                                          
	                  Kappa : 0.6241          
	 Mcnemar's Test P-Value : < 2.2e-16       

	Statistics by Class:

	                     Class: archea Class: bact Class: euk Class: viruses
	Sensitivity                 0.6401      0.8528     0.7750         0.5684
	Specificity                 0.9415      0.8728     0.8702         0.9384
	Pos Pred Value              0.7528      0.7035     0.7268         0.7137
	Neg Pred Value              0.9039      0.9437     0.8967         0.8894
	Prevalence                  0.2177      0.2614     0.3081         0.2128
	Detection Rate              0.1393      0.2229     0.2388         0.1210
	Detection Prevalence        0.1851      0.3168     0.3286         0.1695
	Balanced Accuracy           0.7908      0.8628     0.8226         0.7534

On the non-reduced version: 
	> confusionMatrix(predicted_classes_full,testClass)
	Confusion Matrix and Statistics

	          Reference
	Prediction archea bact euk viruses
	   archea     388   33  49      47
	   bact        83  636  76      91
	   euk        105   44 688     109
	   viruses     52   41  76     367

	Overall Statistics
	                                          
	               Accuracy : 0.7206          
	                 95% CI : (0.7039, 0.7369)
	    No Information Rate : 0.3081          
	    P-Value [Acc > NIR] : < 2.2e-16       
	                                          
	                  Kappa : 0.6225          
	 Mcnemar's Test P-Value : 2.99e-14        

	Statistics by Class:

	                     Class: archea Class: bact Class: euk Class: viruses
	Sensitivity                 0.6178      0.8435     0.7739         0.5977
	Specificity                 0.9428      0.8827     0.8707         0.9256
	Pos Pred Value              0.7505      0.7178     0.7273         0.6847
	Neg Pred Value              0.8986      0.9410     0.8963         0.8948
	Prevalence                  0.2177      0.2614     0.3081         0.2128
	Detection Rate              0.1345      0.2205     0.2385         0.1272
	Detection Prevalence        0.1792      0.3071     0.3279         0.1858
	Balanced Accuracy           0.7803      0.8631     0.8223         0.7617



# Trying a GC-wise PCA transformation 

We bin the contigs by their GC content, for each bin, we do a PCA with 32 components, we check whether the transformation still capture GC as a first component 

Nothing convincing


# Tring a kNN over the PCA boxCox transformed 5-mers data

Following this blog post (at the bottom)

	http://tgmstat.wordpress.com/2013/11/28/computing-and-visualizing-pca-in-r/

We do 80% training 

We get reasonably good performances when we keep the first 50 dimensions of the PCA :
	> knn_performances
	k-Nearest Neighbors 

	13562 samples
	   50 predictors
	    4 classes: 'archea', 'bact', 'euk', 'viruses' 

	No pre-processing
	Resampling: Cross-Validated (2 fold) 

	Summary of sample sizes: 6781, 6781 

	Resampling results across tuning parameters:

	  k   Accuracy  Kappa  Accuracy SD  Kappa SD
	  11  0.731     0.634  0.00647      0.00906 
	  13  0.734     0.637  0.00542      0.00766 
	  15  0.734     0.637  0.00271      0.00411 
	  17  0.729     0.631  0.00313      0.00485 
	  19  0.73      0.632  0.00407      0.00593 

	Accuracy was used to select the optimal model using  the largest value.
	The final value used for the model was k = 15. 
	> 
	> predicted_classes=predict(knn_performances,testExpr)
	> confusionMatrix(predicted_classes,testClass)
	Confusion Matrix and Statistics

	          Reference
	Prediction archea bact euk viruses
	   archea     479   39  42      47
	   bact       100  844  73     135
	   euk         88   47 850     131
	   viruses     35   30  69     379

	Overall Statistics
	                                          
	               Accuracy : 0.7532          
	                 95% CI : (0.7384, 0.7677)
	    No Information Rate : 0.3052          
	    P-Value [Acc > NIR] : < 2.2e-16       
	                                          
	                  Kappa : 0.6643          
	 Mcnemar's Test P-Value : < 2.2e-16       

	Statistics by Class:

	                     Class: archea Class: bact Class: euk Class: viruses
	Sensitivity                 0.6823      0.8792     0.8221         0.5477
	Specificity                 0.9523      0.8731     0.8870         0.9503
	Pos Pred Value              0.7891      0.7326     0.7616         0.7388
	Neg Pred Value              0.9198      0.9481     0.9190         0.8911
	Prevalence                  0.2072      0.2834     0.3052         0.2043
	Detection Rate              0.1414      0.2491     0.2509         0.1119
	Detection Prevalence        0.1792      0.3400     0.3294         0.1514
	Balanced Accuracy           0.8173      0.8762     0.8545         0.7490



On the topmost 100 PC 

	> knn_performances
	k-Nearest Neighbors 

	13562 samples
	  100 predictors
	    4 classes: 'archea', 'bact', 'euk', 'viruses' 

	No pre-processing
	Resampling: Cross-Validated (2 fold) 

	Summary of sample sizes: 6781, 6781 

	Resampling results across tuning parameters:

	  k   Accuracy  Kappa  Accuracy SD  Kappa SD
	  11  0.724     0.624  0.00209      0.00266 
	  13  0.723     0.622  0.00355      0.00438 
	  15  0.718     0.616  0.00188      0.00242 
	  17  0.721     0.619  0.00156      0.00182 
	  19  0.717     0.613  0.00104      0.0011  

	 > confusionMatrix(predicted_classes,testClass)
	Confusion Matrix and Statistics

	          Reference
	Prediction archea bact euk viruses
	   archea     471   35  47      54
	   bact       107  850  74     131
	   euk         89   36 833     121
	   viruses     35   39  80     386

	Overall Statistics
	                                          
	               Accuracy : 0.7497          
	                 95% CI : (0.7348, 0.7642)
	    No Information Rate : 0.3052          
	    P-Value [Acc > NIR] : < 2.2e-16       
	                                          
	                  Kappa : 0.6599          
	 Mcnemar's Test P-Value : < 2.2e-16       

	Statistics by Class:

	                     Class: archea Class: bact Class: euk Class: viruses
	Sensitivity                 0.6709      0.8854     0.8056         0.5578
	Specificity                 0.9494      0.8715     0.8955         0.9429
	Pos Pred Value              0.7759      0.7315     0.7720         0.7148
	Neg Pred Value              0.9169      0.9506     0.9129         0.8926
	Prevalence                  0.2072      0.2834     0.3052         0.2043
	Detection Rate              0.1390      0.2509     0.2459         0.1139
	Detection Prevalence        0.1792      0.3430     0.3185         0.1594
	Balanced Accuracy           0.8102      0.8785     0.8506         0.75


If we drop the 1st PC : 

	> knn_performances_noPC1
	k-Nearest Neighbors 

	13562 samples
	   99 predictors
	    4 classes: 'archea', 'bact', 'euk', 'viruses' 

	No pre-processing
	Resampling: Cross-Validated (2 fold) 

	Summary of sample sizes: 6781, 6781 

	Resampling results across tuning parameters:

	  k   Accuracy  Kappa  Accuracy SD  Kappa SD
	  11  0.71      0.606  0.00292      0.00453 
	  13  0.708     0.602  0.00156      0.00279 
	  15  0.714     0.611  0.00615      0.0088  
	  17  0.713     0.609  0.00907      0.013   
	  19  0.712     0.608  0.0124       0.0176  

	Accuracy was used to select the optimal model using  the largest value.
	The final value used for the model was k = 15. 
	> confusionMatrix(predicted_classes_no_PC1,testClass)
	Confusion Matrix and Statistics

	          Reference
	Prediction archea bact euk viruses
	   archea     483   43  59      42
	   bact        98  846  89     138
	   euk         90   45 809     127
	   viruses     31   26  77     385

	Overall Statistics
	                                          
	               Accuracy : 0.7447          
	                 95% CI : (0.7296, 0.7593)
	    No Information Rate : 0.3052          
	    P-Value [Acc > NIR] : < 2.2e-16       
	                                          
	                  Kappa : 0.6531          
	 Mcnemar's Test P-Value : < 2.2e-16       

	Statistics by Class:

	                     Class: archea Class: bact Class: euk Class: viruses
	Sensitivity                 0.6880      0.8812     0.7824         0.5564
	Specificity                 0.9464      0.8661     0.8887         0.9503
	Pos Pred Value              0.7703      0.7225     0.7554         0.7418
	Neg Pred Value              0.9207      0.9486     0.9029         0.8930
	Prevalence                  0.2072      0.2834     0.3052         0.2043
	Detection Rate              0.1426      0.2497     0.2388         0.1136
	Detection Prevalence        0.1851      0.3456     0.3161         0.1532
	Balanced Accuracy           0.8172      0.8737     0.8355         0.7533



The SVM-RBF-BoxCox-100-PCA works really well: 

	> svm_performances
	Support Vector Machines with Radial Basis Function Kernel 

	13562 samples
	  100 predictors
	    4 classes: 'archea', 'bact', 'euk', 'viruses' 

	No pre-processing
	Resampling: Cross-Validated (2 fold) 

	Summary of sample sizes: 6781, 6781 

	Resampling results across tuning parameters:

	  sigma    C    Accuracy  Kappa  Accuracy SD  Kappa SD
	  0.00909  0.1  0.645     0.508  0.00271      0.00425 
	  0.00909  1    0.731     0.635  0.0024       0.00404 
	  0.00909  10   0.74      0.649  0.00761      0.00964 
	  0.00909  20   0.738     0.647  0.00574      0.00706 
	  0.01     0.1  0.641     0.502  0.00615      0.00913 
	  0.01     1    0.732     0.637  0.00334      0.00531 
	  0.01     10   0.74      0.648  0.00688      0.00857 
	  0.01     20   0.739     0.648  0.00448      0.00532 
	  0.0111   0.1  0.632     0.489  0.00448      0.0068  
	  0.0111   1    0.734     0.639  0.00334      0.00525 
	  0.0111   10   0.741     0.65   0.00532      0.00645 
	  0.0111   20   0.741     0.65   0.00553      0.00674 
	  0.0125   0.1  0.619     0.47   0.00355      0.00546 
	  0.0125   1    0.735     0.64   0.00334      0.0053  
	  0.0125   10   0.742     0.651  0.00417      0.00496 
	  0.0125   20   0.743     0.652  0.00417      0.00496 

	Accuracy was used to select the optimal model using  the largest value.
	The final values used for the model were sigma = 0.0125 and C = 20. 
	> confusionMatrix(predicted_classes,testClass)
	Confusion Matrix and Statistics

	          Reference
	Prediction archea bact euk viruses
	   archea     536   50  72      54
	   bact        49  776  47      93
	   euk         78   64 835     107
	   viruses     39   70  80     438

	Overall Statistics
	                                          
	               Accuracy : 0.763           
	                 95% CI : (0.7483, 0.7772)
	    No Information Rate : 0.3052          
	    P-Value [Acc > NIR] : < 2e-16         
	                                          
	                  Kappa : 0.6798          
	 Mcnemar's Test P-Value : 0.05329         

	Statistics by Class:

	                     Class: archea Class: bact Class: euk Class: viruses
	Sensitivity                 0.7635      0.8083     0.8075         0.6329
	Specificity                 0.9345      0.9222     0.8942         0.9299
	Pos Pred Value              0.7528      0.8041     0.7703         0.6986
	Neg Pred Value              0.9380      0.9241     0.9136         0.9080
	Prevalence                  0.2072      0.2834     0.3052         0.2043
	Detection Rate              0.1582      0.2290     0.2465         0.1293
	Detection Prevalence        0.2102      0.2848     0.3200         0.1851
	Balanced Accuracy           0.8490      0.8652     0.8509         0.7814
	> 

We save the model as "svm-rbf-pca-90.Rdat"