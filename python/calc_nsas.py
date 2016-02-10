import math as m
import sys

k1 = float(sys.argv[1])
a1 = m.exp(float(sys.argv[2]))/1e10
k2 = float(sys.argv[3])
a2 = m.exp(float(sys.argv[4]))/1e10

print "({:.4g}, {:.4g})".format(k1, a1)
print "({:.4g}, {:.4g})".format(k2, a2)

ns = 1 + m.log(a2/a1)/m.log(k2/k1)
As = a1*(0.05/k1)**(ns-1)

print "ns:", ns
print "As:", As
