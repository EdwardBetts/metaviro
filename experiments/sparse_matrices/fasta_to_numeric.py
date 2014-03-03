import sys
from string import maketrans


intab = "ACGTacgt"
outtab = "01230123"
transtab = maketrans(intab, outtab)

inputfile  = sys.argv[1]
outputfile = sys.argv[2]


with open(inputfile, "r") as inputfh:
	with open(outputfile, "w") as outputfh:
		for line in inputfh.readlines():
			if line[0] != ">":
				outputfh.write(line.translate(transtab))
			else:
				outputfh.write(line)
				
