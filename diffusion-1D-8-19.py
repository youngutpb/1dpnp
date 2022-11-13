#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 23:55:21 2022

@author: douglas
"""
import numpy as np
import matplotlib.pyplot as plt







Nx = 10
dx = 1.0/Nx

D = 0.01

## Stability Criterion:  dt <= dx^2/(2*D)
## https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture16%20--%20Numerical%20methods%20for%20diffusion%20models.html

Nt = 20
#dt = 0.25
dt_min = dx**2/(2*D)
print(dt_min)
dt = dt_min
Ttotal = Nt*dt

x = np.linspace(0,1.0, Nx+1)
print(x)

c = np.zeros((Nt+1, Nx+1))

for ix in np.arange(0, Nx):
    c[0,ix] = 4.0*x[ix]*(1.0-x[ix]) 

#c[0,int(Nx/2-1)] = 1.0
# N_half = int(Nx/2)
# print(N_half)
#c[0,N_half] = 1.0
#c[0,int(.5/dx)] = 1.0

plt.plot(x[:],c[0,:])



for it in range(1, Nt) :
    #print("it = ", it)
    for ix in np.arange(1, Nx) :
        #print(ix)
        if (ix == 1):
            c[it,ix] = c[it-1,ix] + (D*(c[it-1,ix+1]-c[it-1,ix])/dx**2)*dt
            c[it,0] = c[it,1]
            #print("left")
        elif (ix == Nx-1) :
            c[it,ix] = c[it-1,ix] - (D*(c[it-1,ix]-c[it-1,ix-1])/dx**2)*dt
            c[it,Nx] = c[it,Nx-1]
            #print("right")
        else :
            c[it,ix] = c[it-1,ix] + (D*(c[it-1,ix+1]-2*c[it-1,ix]+c[it-1,ix-1])/dx**2)*dt
            #print("mid")
    plt.plot(x[:],c[it,:])
plt.show()