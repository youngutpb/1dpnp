#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 23:09:13 2022

@author: douglas
"""
import numpy as np
import matplotlib.pyplot as plt

Nx = 40
gp_int = Nx -1

#rho_o = 1.0/8.854e-12
rho_o = 1100
#epsilon_o = 1.0 # Electrostatic Constant
#epsilon_o = 8.854E-12 # Electrostatic Constant
Lx = 0.01
x = np.linspace(0,Lx, Nx+1)
print(x)

dx = Lx/(Nx)

matrix = np.zeros((Nx+1, Nx+1))
inverse_matrix = np.zeros((Nx+1, Nx+1))

rho = np.zeros((Nx+1))
#rho[0] = 1.0

phi = np.zeros((Nx+1))

phi_left = 0.0
phi_right = 1.0



rho[0] = rho_o*dx**2 - phi_left
rho[Nx-1] = rho_o*dx**2 - phi_right 

for i in np.arange(1, Nx -1):
    #rho[i] = 4.0*x[i]*(1.0-x[i])
    rho[i] = rho_o*dx**2
    

for j in np.arange(0,Nx -1):
    for i in np.arange(0,Nx -1):
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
    
for j in np.arange(0, gp_int):
    for i in np.arange(0, gp_int):
        if (i >= j):
            element = -(j+1)*(gp_int-i)/(gp_int+1)
        elif (i < j):
            element = -(i+1)*(gp_int-j)/(gp_int+1)
        inverse_matrix[i,j] = element
    
print(inverse_matrix)
print(rho)

# phi = np.linalg.solve(matrix, rho)
# print(phi)

#res = matrix @ inverse_matrix

phi = inverse_matrix @ rho
print(phi)

phi_plot = np.zeros((Nx + 1))
phi_plot[0]= phi_left
phi_plot[Nx] = phi_right

for i in np.arange(1, Nx):
    phi_plot[i] = phi[i-1]
    

phi_plot_analytic = np.zeros((Nx+1))

C_analytic = (phi_right - phi_left - 0.5*rho_o*Lx**2)/Lx
for i in np.arange(0, Nx+1):
    x[i] = i*dx
#    phi_plot_analytic[i] = 0.5*rho[i]*x[i]**2 + C_analytic*x[i] + phi_left
    phi_plot_analytic[i] = (0.5*rho[i]/dx**2)*(x[i]**2 - Lx*x[i])

plt.xlim(0, Lx)
plt.plot(x,phi_plot, label = "numeric")
plt.plot(x,phi_plot_analytic, label = "analytic", linestyle = "--")
plt.legend()
plt.show

    