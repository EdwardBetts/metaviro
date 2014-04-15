#!/bin/bash

OLDIFS=$IFS
IFS=$'\n'
for line in $(cat $1 | grep '^>\|Query=')
do
	
	if [ "$prevline" != "" ]
	then
		if [ "$(echo $prevline | head -c 6)" == "Query=" ]
		then
			if [ "$(echo $line | head -c 1)" == ">" ]
			then
				echo $prevline
				#echo $line
			fi
		fi
	fi
	prevline=$line
done
IFS=$OLDIFS
