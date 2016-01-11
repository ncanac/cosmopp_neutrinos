#from matplotlib import pyplot as plt
import numpy as np
import sys

root = sys.argv[1] # directory + root for data files
par1 = sys.argv[2] # y-axis parameter
par2 = sys.argv[3] # x-axis parameter

params_fname = root + ".paramnames"
posterior_fname = root + "posterior.txt"

paramnames = []
f = open(params_fname, 'r')
for line in f:
    paramnames.append(line.split()[0])
f.close()

par1_idx = paramnames.index(par1)
par2_idx = paramnames.index(par2)
lnlike_idx = len(paramnames)

posterior = []
f = open(posterior_fname, 'r')
for line in f:
    vals = line.split()
    posterior.append([float(vals[par1_idx]), float(vals[par2_idx]), float(vals[lnlike_idx])])
f.close()

posterior = np.array(posterior)
