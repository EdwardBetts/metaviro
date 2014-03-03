import numpy as np
import scipy as sp
import random
import pylab as pl
from matplotlib.colors import ListedColormap
from sklearn import neighbors, datasets

from sklearn.metrics import confusion_matrix, f1_score, SCORERS
from kappa import kappa



n_neighbors = 15

# import some data to play with
iris = datasets.load_iris()
X = iris.data[:, :2]  # we only take the first two features. We could
                      # avoid this ugly slicing by using a two-dim dataset
y = iris.target

h = .02  # step size in the mesh

# Create color maps
cmap_light = ListedColormap(['#FFAAAA', '#AAFFAA', '#AAAAFF'])
cmap_bold = ListedColormap(['#FF0000', '#00FF00', '#0000FF'])

for weights in ['uniform', 'distance']:
    # we create an instance of Neighbours Classifier and fit the data.
    clf = neighbors.KNeighborsClassifier(n_neighbors, weights=weights)
    clf.fit(X, y)

    # Plot the decision boundary. For that, we will assign a color to each
    # point in the mesh [x_min, m_max]x[y_min, y_max].
    x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
    y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                         np.arange(y_min, y_max, h))
    Z = clf.predict(np.c_[xx.ravel(), yy.ravel()])

    # Put the result into a color plot
    Z = Z.reshape(xx.shape)
    pl.figure()
    pl.pcolormesh(xx, yy, Z, cmap=cmap_light)

    # Plot also the training points
    pl.scatter(X[:, 0], X[:, 1], c=y, cmap=cmap_bold)
    pl.xlim(xx.min(), xx.max())
    pl.ylim(yy.min(), yy.max())
    pl.title("3-Class classification (k = %i, weights = '%s')"
             % (n_neighbors, weights))

pl.show()





### Classes of training data, confusion matrix, kappa

training_classes=clf.predict(X)
results=np.matrix((training_classes,y)).transpose()


def my_confusion_matrix(two_col_matrix):
	""""Classes are assumed to be int"""
	dim=max(two_col_matrix[:,0])[0,0]+1
	matrix=np.zeros((dim,dim),np.dtype('i4'))
	for m in two_col_matrix:
		matrix[m[0,0],m[0,1]]+=1

	return matrix


print(my_confusion_matrix(results))
cm=confusion_matrix(results[:,0],results[:,1])
print(cm)
pl.matshow(cm)
pl.title('Confusion matrix')
pl.colorbar()
pl.ylabel('True label')
pl.xlabel('Predicted label')
pl.show()
print(kappa(results[:,0],results[:,1]))



### Distance or uniform? LA says distance, based on "paire de plots"
uniform_model=  neighbors.KNeighborsClassifier(n_neighbors, weights="uniform")
uniform_model.fit(X, y)
uniform_model_training_classes=uniform_model.predict(X)
print(confusion_matrix(uniform_model_training_classes,y))
print(kappa(y,uniform_model_training_classes))

distance_model=  neighbors.KNeighborsClassifier(n_neighbors, weights="distance")
distance_model.fit(X, y)
distance_model_training_classes=distance_model.predict(X)
print(confusion_matrix(y,distance_model_training_classes))
print(kappa(y,distance_model_training_classes))



### Which k? 

for k in range(1,16,1):
	distance_model=  neighbors.KNeighborsClassifier(k, weights="distance")
	distance_model.fit(X, y)
	distance_model_training_classes=distance_model.predict(X)
	# print(confusion_matrix(distance_model_training_classes,y))
	print k,kappa(y,distance_model_training_classes)


### Why for k=1 kappa is not 1 ?

distance_model_k1=  neighbors.KNeighborsClassifier(1, weights="distance")
distance_model_k1.fit(X, y)
distance_model_k1_training_classes=distance_model_k1.predict(X)
# print(confusion_matrix(distance_model_training_classes,y))
print k,kappa(y,distance_model_training_classes)
disagrees=np.nonzero(distance_model_k1_training_classes!=y)
print y[disagrees],distance_model_k1.predict(X[disagrees])


x_min, x_max = X[:, 0].min() - 1, X[:, 0].max() + 1
y_min, y_max = X[:, 1].min() - 1, X[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, h),
                     np.arange(y_min, y_max, h))
Z = distance_model_k1.predict(np.c_[xx.ravel(), yy.ravel()])

# Put the result into a color plot
Z = Z.reshape(xx.shape)


pl.subplot(1,3,1)
pl.pcolormesh(xx, yy, Z, cmap=cmap_light)
pl.scatter(X[:, 0], X[:, 1], c=y,cmap=cmap_bold)
 
pl.subplot(1,3,2)
pl.pcolormesh(xx, yy, Z, cmap=cmap_light)
pl.scatter(X[disagrees, 0], X[disagrees, 1], c=y[disagrees],cmap=cmap_bold,marker="x")


pl.subplot(1,3,3)
pl.pcolormesh(xx, yy, Z, cmap=cmap_light)
pl.scatter(X[disagrees, 0], X[disagrees, 1], c=distance_model_k1.predict(X[disagrees]),cmap=cmap_bold,marker="x")
pl.show()



### Hyp: Several observations are equals, explaining why kappa != 1 for 1-NN 
equals=np.nonzero(np.bitwise_and(X[:,0]==6.4, X[:,1]==3.2))
print np.column_stack((X[equals],y[equals]))



### Take all four features
X = iris.data[:, :4]  # we only take the first two features. We could
                      # avoid this ugly slicing by using a two-dim dataset

distance_model_k1=  neighbors.KNeighborsClassifier(1, weights="distance")
distance_model_k1.fit(X, y)
distance_model_k1_training_classes=distance_model_k1.predict(X)
# print(confusion_matrix(distance_model_training_classes,y))
print 1,kappa(y,distance_model_k1_training_classes)
disagrees=np.nonzero(distance_model_k1_training_classes!=y)

if disagrees[0].size>0:
	print y[disagrees],distance_model_k1.predict(X[disagrees])

## We can't plot the decision boundaries in the 4D space, but we can project
### Not enough memory, skip


### We verify that every points are unique with the four features
### We create a distance matrix 
pairwise_dist=sp.spatial.distance.pdist(X)
print min(pairwise_dist)
equal_features=np.nonzero(pairwise_dist==0)

# Coordinates of equal vectors 
pairwise_dist_squares=sp.spatial.distance.squareform(pairwise_dist)
equal_vectors=np.nonzero(pairwise_dist_squares==0)
for r in np.column_stack(equal_vectors):
	if r[0]!=r[1]:
		print r[0],r[1]

# We print the coord of some equal vec
print X[142,:],y[142]
print X[101,:],y[101]




distance_model_k15=  neighbors.KNeighborsClassifier(15, weights="distance")
distance_model_k15.fit(X, y)
distance_model_k15_training_classes=distance_model_k15.predict(X)
print(confusion_matrix(distance_model_k15_training_classes,y))
print 15,kappa(y,distance_model_k15_training_classes)

## For k==15, 4 features, we have kappa=1

### What happens if we take a random subset for training and test on the rest? 


train_index=random.sample(range(1,len(X)),20)
test_index=list(set(range(1,len(X))).difference(train_index))

x_train=X[train_index]
y_train=y[train_index]
x_test=X[test_index]
y_test=y[test_index]

distance_model_k15_test_train=  neighbors.KNeighborsClassifier(15, weights="distance")
distance_model_k15_test_train.fit(x_train, y_train)
distance_model_k15_test_train_training_classes=distance_model_k15_test_train.predict(x_test)
print(confusion_matrix(distance_model_k15_test_train_training_classes,y_test))
print 15,kappa(y_test,distance_model_k15_test_train_training_classes)


#### Plot the distribution of kappa as a function of number of neighbors (k) and the size of the training set 
