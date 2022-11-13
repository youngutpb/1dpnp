#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 23:09:13 2022

@author: douglas
"""
import numpy as np
import matplotlib.pyplot as plt

gridpoints = 20

x = np.linspace(0,1.0, gridpoints)
print(x)

matrix = np.zeros((gridpoints, gridpoints))
rho = np.zeros((gridpoints))
rho[0] = 1.0

phi = np.zeros((gridpoints))

for j in np.arange(0,gridpoints ):
    for i in np.arange(0,gridpoints):
        diff = i - j
        if (diff == -1):
            element = -1.0
        elif (diff == 1):
            element = 1.0
        elif (diff == 0) :
            element = 2.0
        else:
            element = 0.0
        matrix[i,j] = element
        
    #rho[j] = 4.0*x[j]*(1.0-x[j]) 
    
print(matrix)
print(rho)

phi = np.linalg.solve(matrix, rho)
print(phi)
print(np.shape(phi))

plt.plot(x,phi)
plt.show

    