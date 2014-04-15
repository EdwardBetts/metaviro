#!/bin/bash

for filename in sample_matrix_*.bdat
do
	python ../knn.py $filename >> tagadatsouintsouin.log
done
