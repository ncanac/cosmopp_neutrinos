from matplotlib import pyplot as plt
import numpy as np
import sys
from matplotlib.patches import Ellipse
from scipy.interpolate import griddata

root = sys.argv[1] # directory + root for data files
xpar = sys.argv[2] # x-axis parameter
#xmin = float(sys.argv[3]) # minimum value of x parameter
#xmax = float(sys.argv[4]) # maximum value of x parameter
ypar = sys.argv[3] # y-axis parameter
#ymin = float(sys.argv[6]) # minimum value of y parameter
#ymax = float(sys.argv[7]) # maximum value of y parameter
ndiv = float(sys.argv[4])

params_fname = root + ".paramnames"
samples_fname = root + "posterior.txt"
constraints_fname = root + "parameter_constraints.txt"

paramnames = []
f = open(params_fname, 'r')
for line in f:
    paramnames.append(line.split()[0])
f.close()

xpar_idx = paramnames.index(xpar)
ypar_idx = paramnames.index(ypar)
lnlike_idx = len(paramnames)+1

samples = []
f = open(samples_fname, 'r')
for line in f:
    vals = line.split()
    samples.append([float(vals[xpar_idx]), float(vals[ypar_idx]), float(vals[lnlike_idx])])
f.close()

samples = np.array(samples)

f = open(constraints_fname, 'r')
for line in f:
    vals = line.split()
    if xpar in vals[0]:
        i = vals[2].find("+-")
        xmean = float(vals[2][:i])
        xsd = float(vals[2][i+2:])
        xmin = max(0.0, xmean - 3.0*xsd)
        xmax = xmean + 3.0*xsd
    elif ypar in vals[0]:
        i = vals[2].find("+-")
        ymean = float(vals[2][:i])
        ysd = float(vals[2][i+2:])
        ymin = max(0.0, ymean - 3.0*ysd)
        ymax = ymean + 3.0*ysd

#ndiv = 40.0
xbin = (xmax - xmin) / ndiv
ybin = (ymax - ymin) / ndiv
x_grid, y_grid = np.meshgrid(np.arange(xmin, xmax, xbin), np.arange(ymin, ymax, ybin))

prob = np.zeros_like(x_grid)
counts = np.zeros_like(x_grid)

for row in samples:
    xval = row[0]
    yval = row[1]
    probval = row[2]
    if xval > xmin and xval < xmax and yval > ymin and yval < ymax:
        xi = int((xval - xmin) / xbin)
        yi = int((yval - ymin) / ybin)
        prob[xi, yi] += probval
        counts[xi, yi] += 1.0

levels = [0.0]
sorted_prob = sorted(prob.flatten(), reverse=True)
cum_prob = np.cumsum(sorted_prob)
#levels.append(sorted_prob[np.where(cum_prob >= 0.9973)[0][0]])
levels.append(sorted_prob[np.where(cum_prob >= 0.9545)[0][0]])
levels.append(sorted_prob[np.where(cum_prob >= 0.6827)[0][0]])
levels.append(sorted_prob[0])

#fig, ax = plt.subplots(figsize=(10, 7.5))
plt.figure(figsize=(10, 7.5))
maxvalue = np.max(prob)
#levels = [0, maxvalue / 100, maxvalue / 30, maxvalue / 10, maxvalue / 3, maxvalue]
cplot = plt.contourf(x_grid, y_grid, prob, levels=levels)#, 50, cmap="RdBu")#, vmin=0, vmax=1)

cbar = plt.colorbar(cplot)
#cbar.set_label("$2 \Delta ln(L)$", size=20)
#cbar.set_ticks([0, 1.0, 3.0, 9.0, 16.0, 25.0, 1e6])

#plt.xlim(xmin - xbin, xmax + xbin)
#plt.ylim(ymin - ybin, ymax + ybin)
plt.xlabel(xpar, size=20)
plt.ylabel(ypar, size=20)
plt.tight_layout()

plt.show()
#plt.savefig("name.png")
