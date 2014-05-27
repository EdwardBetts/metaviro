import sys, re
from string import maketrans
from Bio import SeqIO
from Bio.Seq import Seq

intab = "ACGTacgt"
outtab = "01230123"
transtab = maketrans(intab, outtab)

inputfile  = sys.argv[1]
outputfile = sys.argv[2]

recordlist = []
with open(inputfile, "r") as inputfh:
	for record in SeqIO.parse(inputfh, "fasta"):
		if re.match('^[ATGC0123]*$', str(record.seq).upper()):
			record.seq = Seq(str(record.seq).translate(transtab))
			recordlist.append(record)
				
with open(outputfile, "w") as outputfh:
	SeqIO.write(recordlist, outputfh, "fasta")
