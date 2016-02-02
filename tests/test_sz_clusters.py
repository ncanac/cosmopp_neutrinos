import numpy as np
from scipy.stats import norm
import sys

filename = sys.argv[1]

mean = 0.78
sd = 0.01

n = norm(loc=mean, scale=sd)

params = []
p_tot = 0
with open(filename, 'r') as f:
    for line in f:
        vals = [float(val) for val in line.split()]
        ombh2 = vals[0]
        omch2 = vals[1]
        h = vals[2]
        p_old = vals[-2]
        sigma8 = vals[-1]
        omm = (ombh2 + omch2)/h**2.
        x = sigma8*(omm/0.27)**0.3
        p_new = p_old*norm.pdf(x)
        p_tot += p_new
        vals[-2] = p_new
        params.append(vals)

params = np.array(params)
params[:,-2] = params[:,-2]/p_tot

assert sum(params[:,-2]) < 1.01 and sum(params[:,-2]) > 0.99

for sample in params:
    s = ""
    for par in sample:
        s += str(par) + " "
    print s
