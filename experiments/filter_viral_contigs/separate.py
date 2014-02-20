
import sys

fastaf = sys.argv[1]
hitsf = sys.argv[2]

with open(hitsf, "r") as hitsfh:
	hits = hitsfh.readlines()

with open(fastaf, "r") as fastafh:
	keep = False
	for line in fastafh.readlines():
		if line[0] == ">":
			if line[1:] not in hits:
				keep = True
				sys.stdout.write(line)
			else:
				keep = False
				#sys.stderr.write("skipping %s" % line)
		else:
			if keep:
				sys.stdout.write(line)


