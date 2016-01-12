#from matplotlib import pyplot as plt
import numpy as np
import sys

root = sys.argv[1] # directory + root for data files
xpar = sys.argv[2] # x-axis parameter
xmin = float(sys.argv[3]) # minimum value of x parameter
xmax = float(sys.argv[4]) # maximum value of x parameter
ypar = sys.argv[5] # y-axis parameter
ymin = float(sys.argv[6]) # minimum value of y parameter
ymax = float(sys.argv[7]) # maximum value of y parameter

params_fname = root + ".paramnames"
posterior_fname = root + "posterior.txt"

paramnames = []
f = open(params_fname, 'r')
for line in f:
    paramnames.append(line.split()[0])
f.close()

xpar_idx = paramnames.index(xpar)
ypar_idx = paramnames.index(ypar)
lnlike_idx = len(paramnames)

posterior = []
f = open(posterior_fname, 'r')
for line in f:
    vals = line.split()
    posterior.append([float(vals[xpar_idx]), float(vals[ypar_idx]), float(vals[lnlike_idx])])
f.close()

posterior = np.array(posterior)

ndiv = 10.0
xspc = (xmax - xmin) / ndiv
yspc = (ymax - ymin) / ndiv
x, y = np.meshgrid(np.arange(xmin, xmax, xspc), np.arange(ymin, ymax, yspc))

likelihood = np.zeros_like(x)
likelihood_counts = np.zeros_like(x)

for row in posterior:
    xval = row[0]
    yval = row[1]
    likeval = row[2]
    xi = int((xval - xmin) / xspc)
    yi = int((yval - ymin) / yspc)
    likelihood[xi, yi] += likeval
    likelihood_counts[xi, yi] += 1.0

likelihood = likelihood / likelihood_counts
min_like = np.min(likelihood)
likelihood = likelihood - min_like

fig, ax = plt.subplots(figsize=(10, 7.5))
contour = ax.contourf(xx, yy, probs, 50, cmap="RdBu", vmin=0, vmax=1)

ax_c = f.colorbar(contour)
ax_c.set_label("$2 \Delta ln(L)$", size=20)
#ax_c.set_ticks([0, 1.0, 3.0, 9.0, 16.0, 25.0, 1e6])

plt.xlim(x1min - x1pad, x1max + x1pad)
plt.ylim(x2min - x2pad, x2max + x2pad)
plt.xlabel("$" + xpar + "$", size=20)
plt.ylabel("$" + ypar + "$", size=20)
plt.tight_layout()

plt.show()
#plt.savefig("log_reg_nl.png")
