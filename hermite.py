#!/usr/bin/env python

import matplotlib.pyplot as plt
import math
import numpy as np
import sys
import argparse

def H(x,n):
    if n == 0 :
        return 1.0
    elif n == 1 :
        return 2.0*x
    elif n == 2 :
        return 4.0*x**2-2.0
    elif n == 3 :
        return 8.0*x**3-12.0*x
    else :
        print "error"
        sys.exit(1)

def hermite(n,x):
    A = math.sqrt(math.sqrt(0.5/math.pi)/(2**n*math.factorial(n)))
    gsi = math.sqrt(0.5)*x
    return A*H(gsi,n)*math.exp(-0.5*gsi**2)

def plot(n):
    # draw positive x
    x = np.array(range(0,100),dtype=float)
    x *= 0.1
    y = [hermite(n,i) for i in x]
    plt.plot(x,y,"r-")

    # draw negative x
    x = np.array(range(0,-100,-1),dtype=float)
    x *= 0.1
    y = [hermite(n,i) for i in x ]
    print x
    print y
    plt.plot(x,y,"r-")

    # plot
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-n","--n",type=int,required=True)
    args = parser.parse_args()
    plot(args.n)
