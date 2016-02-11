import sys

root = sys.argv[1]

with open(root + "posterior_sigma8.txt", 'r') as infile:
    with open(root + "sigma8_posterior.txt", 'w') as outfile:
        for line in infile:
            vals = line.split()
            s = vals[-1]
            for v in vals[0:-1]:
                s += " " + v
            outfile.write(s + "\n")
