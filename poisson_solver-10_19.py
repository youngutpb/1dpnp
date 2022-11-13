#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 23:09:13 2022

@author: douglas
"""
import numpy as np
import matplotlib.pyplot as plt

Nx = 200

epsilon_o = 8.854E-12 # Electrostatic Constant
cn = 5.0E12 # concentration/number density in m^-3
q_pro = 1.602e-19 # charge on a proton
rho_o = cn*q_pro
rho_epo = rho_o/epsilon_o
#epsilon_o = 1.0 # Electrostatic Constant

Lx = 0.01
x = np.linspace(0,Lx, Nx+1)
print(x)

dx = Lx/(Nx)

matrix = np.zeros((Nx+1, Nx+1)) # Poisson Matrix for check
inverse= np.zeros((Nx+1, Nx+1)) # Inverse of Poission Matrix

rho = np.zeros((Nx+1))
#rho[0] = 1.0

phi = np.zeros((Nx+1))

phi_left = 0.0
phi_right = 1.0



rho[0] = rho_epo*dx**2 - phi_left
rho[Nx] = rho_o*dx**2 - phi_right 

for i in np.arange(1, Nx):
    #rho[i] = 4.0*x[i]*(1.0-x[i])
    rho[i] = rho_epo*dx**2
    

for j in np.arange(0,Nx +1):
    for i in np.arange(0,Nx +1):
        diff = i - j
        if (diff == -1):
            element = 1.0
        elif (diff == 1):
            element = 1.0
        elif (diff == 0) :
            element = -2.0
        else:
            element = 0.0
        matrix[i,j] = element
        
    #rho[j] = 4.0*x[j]*(1.0-x[j]) 
    
for j in np.arange(0, Nx+1):
    for i in np.arange(0, Nx+1):
        if (i >= j):
            element = -(j+1)*(Nx+1-i)/(Nx+1+1)
        elif (i < j):
            element = -(i+1)*(Nx+1-j)/(Nx+1+1)
        inverse[i,j] = element
    
print(inverse)
print(rho)

# phi = np.linalg.solve(matrix, rho)
# print(phi)

res = matrix @ inverse
print(res)

phi = inverse @ rho
print(phi)

phi_plot = np.zeros((Nx + 1))
# phi_plot[0]= phi_left
# phi_plot[Nx] = phi_right

for i in np.arange(1, Nx):
    phi_plot[i] = phi[i]
    

phi_plot_analytic = np.zeros((Nx+1))
C_analytic = phi_right/Lx - 0.5*rho_epo*Lx - phi_left/Lx
for i in np.arange(0, Nx+1):
    x[i] = i*dx
    #phi_plot_analytic[i] = 0.5*rho[i]*x[i]**2 + C_analytic*x[i] + phi_left # Analytic solution for rho = 0
    #phi_plot_analytic[i] = (0.5*rho[i])*(x[i]**2 - Lx*x[i]) # Analytic solution for rho = constant  <=  NEED TO FIX THIS FOR ARB PHI_RIGHT & PHI_LEFT
    phi_plot_analytic[i] = 0.5*rho_epo*x[i]**2 + C_analytic*x[i] + phi_left

plt.xlim(0, Lx)
plt.plot(x,phi, label = "numeric")
plt.plot(x,phi_plot_analytic, label = "analytic", linestyle = "--")
plt.legend()
plt.show

    