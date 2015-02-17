# how to detect wildcards pumping up

#{sample}.txt

#1.txt
#11.txt
#111.txt
#1111.txt


#1.a.txt
#1.a.1.txt
#1.a.1.a.txt


#check wildcard value for periodicity?
#how to implement this?
#have a look at regular languages. Is there an automaton that detects this?
import re
repeat = re.compile("^((?P<value>.{1,50})(?P=value){9,100})$")
repeat = re.compile("^((?P<value>.+)(?P=value){9,100})$")

v = "1a1a1a1a1a1a1a1a1a1a"

print(repeat.match(v).groups())
