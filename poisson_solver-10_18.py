#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 19 23:09:13 2022

@author: douglas
"""
import numpy as np
import matplotlib.pyplot as plt

Nx = 100
gp_int = Nx +1

#rho_o = 1.0/8.854e-12
#rho_o = 0.0
rho_o = 100000
Lx = 0.01
x = np.linspace(0,Lx, Nx+1)
#print(x)

dx = Lx/(Nx)

matrix = np.zeros((Nx+1, Nx+1))
inverse_matrix = np.zeros((gp_int, gp_int))

rho = np.zeros((Nx+1))
rho_num = np.zeros((Nx+1))
#rho[0] = 1.0

phi = np.zeros((Nx+1))

phi_left = 2.0
phi_right = 1.0

rho[0] = rho_o
rho[Nx-1] = rho_o

rho_num[0] = rho_o*dx**2 - phi_left
rho_num[Nx] = rho_o*dx**2 - phi_right 

for i in np.arange(1, Nx):
    #rho[i] = 4.0*x[i]*(1.0-x[i])
    rho[i] = rho_o
    rho_num[i] = rho_o*dx**2
    
#=========== Calculates the Poison Matrix ==================
for j in np.arange(0,Nx+1):
    for i in np.arange(0,Nx+1):
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
print("Poisson Matrix =", matrix)    


#=========== Calculates the Inverse of the Poison Matrix ==================
### Reference:  https://www.scirp.org/journal/paperinformation.aspx?paperid=50041
for j in np.arange(0, gp_int):
    for i in np.arange(0, gp_int):
        if (i >= j):
            element = -(j+1)*(gp_int-(i))/(gp_int+1)
        elif (i < j):
            element = -(i+1)*(gp_int-(j))/(gp_int+1)
        inverse_matrix[i,j] = element
    
print("Inverse Matrix =", inverse_matrix)
print("rho_num =", rho_num)

# phi = np.linalg.solve(matrix, rho)
# print(phi)

#res = matrix @ inverse_matrix


phi = inverse_matrix @ rho_num

phi[0] = phi_left
phi[Nx] = phi_right
print("phi =", phi)

phi_plot = np.zeros((Nx + 1))
phi_plot[0]= phi_left
phi_plot[Nx] = phi_right

for i in np.arange(1, Nx+1):
    phi_plot[i] = phi[i-1]
    

phi_plot_analytic = np.zeros((Nx+1))

#C_analytic = (phi_right - phi_left - 0.5*rho_o*Lx**2)/Lx
C_analytic = phi_right/Lx - 0.5*rho_o*Lx - phi_left/Lx
for i in np.arange(0, Nx):
    x[i] = i*dx
    #phi_plot_analytic[i] = 0.5*rho[i]*x[i]**2 + C_analytic*x[i] + phi_left # Analytic solution for rho = 0
    #phi_plot_analytic[i] = (0.5*rho[i])*(x[i]**2 - Lx*x[i]) # Analytic solution for rho = constant  <=  NEED TO FIX THIS FOR ARB PHI_RIGHT & PHI_LEFT
    phi_plot_analytic[i] = 0.5*rho[i]*x[i]**2 + C_analytic*x[i] + phi_left
    
plt.xlim(0, Lx)
plt.plot(x,phi_plot, label = "numeric")
plt.plot(x,phi_plot_analytic, label = "analytic", linestyle = "--")
plt.legend()
plt.show

    