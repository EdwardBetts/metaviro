
import sys

prevline = ''
with open(sys.argv[1], 'r') as blastout:
	for line in blastout.readlines():
		if line[:6] != "Query=" and line[0] != ">":
			continue
		if prevline != '' and prevline[:6] == "Query=" and line[0] != ">":
			sys.stdout.write(prevline[7:])
		prevline = line
if prevline[0] != ">":
	sys.stdout.write(prevline[7:])


