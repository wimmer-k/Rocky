#!/usr/bin/env python2.7
import sys
import numpy as np
import matplotlib.pyplot as plt
#import os

def potential(x):
    return 0

ds = 0.1 #step size
ds2 = ds*ds #step size
n = 10   #n steps
phi = []
psi = []

def numerovstep(x):
    phi[2] = ds2*f*psi[1]+2.*phi[1]-phi[0]
    phi[0] = phi[1]
    phi[1] = phi[2]
    f = 2.*(potential(x)-e)
    psi[1] = phi[1]/(1.-ds2/12*f)

def main(argv):
    print "hallo"


def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

if __name__ == "__main__":
    sys.exit(main(sys.argv))
         
