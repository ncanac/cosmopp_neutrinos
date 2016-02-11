import numpy as np
from matplotlib import pyplot as plt
import sys
import os

# Command line arguments:
# 1 = filename root
root = sys.argv[1]

d = root.rfind("/")
if d >= 0:
    directory = root[:d+1]
    fname = root[d+1:]
else:
    directory = "./"
    fname = root

# Read in the file names of the chains
chains_list = []
chains_list = [directory + f for f in os.listdir(directory) if (fname in f and ".txt" in f)]

# Read in the parameter names from root.paramnames
param_names = []
f = open(root + ".paramnames", 'r')
for line in f:
    param_names.append(line.split()[0])
f.close()

params = [[] for i in range(len(param_names))]
# Loop through each chain
for chain in chains_list:
    f = open(chain, 'r')
    # Loop through each line in file
    for line in f:
        vals = line.split()
        # Loop through each parameter
        for i,par in enumerate(params):
            # Append value of each parameter, ignoring first two
            params[i].append(float(vals[i+2]))
    f.close()

# Plot histogram for each parameter using 50 bins
for i,par in enumerate(params):
    plt.figure()
    plt.title(param_names[i])
    plt.hist(par, 50, histtype='stepfilled')
plt.show()
