"""
Calculates the new sample probability values from root.txt
for a chain file using a prior on Omega_b, sigma_8, and h
and outputs values to stdout.
"""

import numpy as np
from scipy.stats import norm
import sys

filename = sys.argv[1]

# Fixed hydrostatic bias
#mean = 0.78
#sd = 0.01
# Bias varied in range [0.7, 1.0]
mean = 0.764
sd = 0.025

n = norm(loc=mean, scale=sd)

params = []

p_tot = 0
p_norm = 0
with open(filename, 'r') as f:
    for line in f:
        vals = [float(val) for val in line.split()]
        p = vals[0]
        p_tot += p
        ombh2 = vals[2]
        omch2 = vals[3]
        h = vals[4]
        sigma8 = vals[-1]
        omm = (ombh2 + omch2)/h**2.
        x = sigma8*(omm/0.27)**0.3
        p *= norm.pdf(x)
        p_norm += p
        vals[0] = p
        params.append(vals)

assert p_tot < 1.01 and p_tot > 0.99

params = np.array(params)
params[:, 0] = params[:, 0]/p_norm

assert sum(params[:, 0]) < 1.01 and sum(params[:, 0]) > 0.99

for sample in params:
    s = ""
    for par in sample:
        s += "    " + str(par)
    print s
