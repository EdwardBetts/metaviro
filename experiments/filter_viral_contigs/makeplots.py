"""
Bar chart demo with pairs of bars grouped for easy comparison.
"""
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from collections import defaultdict

csvfile = os.path.splitext(sys.argv[1])[0]

### DATA ###

colors = ['c', 'm', 'y', 'b', 'r', 'g']

print "%s.csv" % csvfile
data=defaultdict(list)
taxa=[]
legend=None
with open("%s.csv" % csvfile) as csvfh:
	legend = csvfh.readline().split()[1:]
	for line in csvfh.readlines():
		line = line.split()
		taxa.append(line[0])
		for i in range(len(legend)):
			try:
				data[legend[i]].append(int(line[i+1]))
			except Exception:
				print "i", i
				print "legend", legend
				print "line", line
				sys.exit(1)

### PLOT ###

n_groups = len(taxa)

fig, ax = plt.subplots()

index = np.arange(n_groups)
bar_width = 0.15

opacity = 0.4
error_config = {'ecolor': '0.3'}
for i in range(len(legend)):
	plt.bar(index + i*bar_width, data[legend[i]], bar_width,
				 alpha=opacity,
				 color=colors[i],
#                     yerr=std_men,
#                     error_kw=error_config,
				 label=legend[i])


plt.xlabel('Taxa')
plt.ylabel('Number of contigs')
plt.title(csvfile)
plt.xticks(index + bar_width, taxa)
plt.legend()

plt.tight_layout()
#    plt.show()
plt.savefig("%s.png" % csvfile)
