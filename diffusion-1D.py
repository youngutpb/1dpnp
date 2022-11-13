#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 23:55:21 2022

@author: douglas
"""
import numpy as np
import matplotlib.pyplot as plt

dx = 0.1
Nt = 10
dt = 0.1
Ttotal = Nt*dt




Nx = int(1/dx)

D = 0.01

t = np.linspace(0, Ttotal, Nt)

x = np.linspace(0,1.0, Nx)

c = np.zeros((Nt, Nx))
# for ix in np.arange(0, Nx):
#     c[0,ix] = 4.0*x[ix]*(1-x[ix]) 

c[0,int(Nx/2-1)] = 1.0
c[0,int(Nx/2)] = 1.0
c[0,int(Nx/2+1)] = 1.0

plt.plot(x[:],c[0,:])



for it in range(1, Nt) :
    print("it = ", it)
    for ix in range(0, Nx) :
        print(ix)
        if (ix == 0):
            c[it,ix] = c[it-1,ix] + (D*(c[it-1,ix+1]-c[it-1,ix])/dx**2)*dt
            print("left")
        elif (ix == Nx-1) :
            c[it,ix] = c[it-1,ix] - (D*(c[it-1,ix]-c[it-1,ix-1])/dx**2)*dt
            print("right")
        else :
            c[it,ix] = c[it-1,ix] + (D*(c[it-1,ix+1]-2*c[it-1,ix]+c[it-1,ix-1])/dx**2)*dt
            print("mid")
    plt.plot(x[:],c[it,:])
plt.show()