from __future__ import division
import sys

fname = sys.argv[1]

badcount = 0
goodcount = 0
with open(fname) as f:
    for line in f:
        count = False
        if "BADHALOFIT" in line:
            badcount += 1
            count = True
        if "GOODHALOFIT" in line:
            goodcount += 1
            count = True
        if count:
            print badcount/(badcount + goodcount)
