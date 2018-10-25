#!/usr/bin/env python
import sys

x = 0
y = 0
lines = 0
for line in sys.stdin:
    items = line.strip().split(",")
    x = x+ float(items[0])
    y = y+ float(items[1])
    lines = lines+1

print "avg: %s - %s" % ( x/lines, y/lines )
