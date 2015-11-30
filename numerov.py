#!/usr/bin/env python2.7
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from array import array

#import os

verbose = 0

#nstep = 402   #n steps
#dx = 5./nstep #step size
dx = 0.01 #step size
dx2 = dx*dx #step size
xmax  = 5
nstep = int(xmax/dx)
start = 0.0

psi = np.zeros(2*nstep)
phi = np.zeros(2*nstep)
x = np.zeros(2*nstep)
k = np.zeros(nstep)
f = np.zeros(nstep)

m=1
hbar =1

def pot(x):
    return 0.5*x*x
#    return -10

def hermite(x,n):
    if n==0:
        return 1
    if n==1:
        return 2*x
    return 2*x*hermite(x,n-1) - 2*(n-1)*hermite(x,n-2)

def exact(nodes):
    for i in range(0,2*nstep):
        phi[i] = math.exp(-x[i]*x[i]/2)*hermite(x[i],nodes)
        
def normalize(wf,norm=1):
    sum =0
    for i in range(len(wf)):
        sum+=wf[i]*wf[i]*dx
    for i in range(len(wf)):
        wf[i]*=math.sqrt(norm/sum)

def generatef(E):
    for i in range(0,nstep):
        x[i] = start+i*dx
        k[i] = 2*m/hbar/hbar*(E-pot(start+i*dx))
        f[i] = 1+dx2/12*k[i]
#        if f[i] == 0:
#            f[i] = 1e-20
            
        if verbose > 1:
            print "%f\t%f\t%f\t%f" % (x[i],pot(start+i*dx),k[i],f[i])

def numerovstep(n):
    psi[n] = ((12-10*f[n-1])*psi[n-1] - f[n-2]*psi[n-2])/f[n]

def main(argv):
    print "hallo"
    nodes = 1

    nsteps = 401
    dx = 5./nsteps
    dx2=dx*dx

    if len(argv)==2:
        nodes = int(argv[1])

    elower = pot(start+(nstep-1)*dx)
    eupper = elower
    for i in range(0,nstep):
        if verbose > 2:
            print "pot(%f+%d*%f) = pot(%f) =  %f" %(start,i,dx,(start+i*dx),pot(start+i*dx))
        if pot(start+i*dx) < elower:
            elower = pot(start+i*dx)
        if pot(start+i*dx) > eupper:
            eupper = pot(start+i*dx)
    
    if elower == eupper:
        elower =0
        eupper =10

    for iter in range(0,1000):
        etrial = 0.5*(elower+eupper)
        if verbose > 2:
            print "%f < E < %f, Etrial = %f" % (elower,eupper,etrial)
        generatef(etrial)
        #if iter ==9:
        #    for i in range(0,nstep):
        #        print "%f\t%f\t%f\t%f" % (x[i],pot(start+i*dx),k[i],f[i])

        #sign of k indicates if the point is in the clasically allowed or forbidden region
        zero_crossings = np.where(np.diff(np.signbit(k)))[0]
        if len(zero_crossings) == 0:
            print "no turning point? unbound state?"
            return
        if verbose >3:
            print zero_crossings[-1]
        if zero_crossings[-1] > nstep-2:
            print "sign change to far, increase x range"
            return

        #initial conditions
        if nodes%2 == 0:
            #even number of nodes, wave function even
            psi[0] = 1
            #assuming f["1"] = f["-1"]
            psi[1] = (12-10*f[0])*psi[0]/2/f[1]
        else:
            #odd wave function
            psi[0] = 0
            psi[1] = dx

        for i in range(2,nstep):
            numerovstep(i)

        zero_crossings = np.where(np.diff(np.signbit(psi)))[0]
        if len(zero_crossings)>0 and zero_crossings[-1] == nsteps:
            zero_crossings = np.delete(zero_crossings, -1)

        if verbose > 1:
            print zero_crossings
            print "iter = %d, E = %f, cross = %d, nodes = %d" %(iter,etrial,len(zero_crossings),int(int(nodes)/2))
        if len(zero_crossings) > int(int(nodes)/2):
            if verbose > 2:
                print "too many nodes found, energy too high"
            eupper = etrial
        else:
            if verbose > 2:
                print "too few nodes found, energy too low"
            elower = etrial

        if (eupper-elower)< 1E-10:
            break

    if verbose >1:
        print "nodes %d, int(int(nodes)/2) %d" %(nodes,int(int(nodes)/2)) 
    print "number of interations: %d, energy: %f" %(iter,etrial)


    for i in range(0,nstep):
        x[nstep+i] = x[i]
        psi[nstep+i] = psi[i]
    for i in range(0,nstep):
        x[i] = -x[2*nstep-i-1]
        if nodes%2==0:
            psi[i] = psi[2*nstep-i-1]
        else:
            psi[i] = -psi[2*nstep-i-1]
        

    file = open('kazz', 'w')
    for i in range(0,nstep*2):
        file.write("%d\t%f\t%f\n" % (i,x[i],psi[i]))
        
    normalize(psi)
    plt.plot(x,psi)
    exact(nodes)
    normalize(phi)
    plt.plot(x,phi)
    plt.show()


def my_range(start, end, step):
    while start <= end:
        yield start
        start += step

if __name__ == "__main__":
    sys.exit(main(sys.argv))
         
