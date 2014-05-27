for filename in *.bdat
do
	python ../make_lmnn_metric.py $filename
done
