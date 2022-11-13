#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 23:55:21 2022

@author: douglas
"""
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import os, subprocess
from matplotlib.ticker import (MultipleLocator,
                               FormatStrFormatter,
                               AutoMinorLocator)

plt.rcParams["font.family"] = "times"
plt.rcParams['text.usetex'] = True

Nx = 40
Lx = 1.0e-2
dx = Lx/Nx
gp_init = Nx - 1


# dx = 1.0e-5
# Lx = 0.2
# Nx = int(Lx/dx)

Nt = 2

# https://www.aqion.de/site/diffusion-coefficients
# https://www.physiologyweb.com/calculators/diffusion_time_calculator.html
Dpos = 1.33e-9 #Na+ m^2/s
Dneg = 2.03e-9 # Cl- m^2/s
c_o = 4.0e13
rho_eo = c_o*1.602e-19/8.854e-12

## Stability Criterion:  dt <= dx^2/(2*D)
## https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture16%20--%20Numerical%20methods%20for%20diffusion%20models.html
#dt = 0.25
dt_min = dx**2/(2*max(Dpos,Dneg))
print(dt_min)
dt = dt_min
Ttotal = Nt*dt

#x = np.linspace(0,Lx, Nx+1)
x = np.zeros((Nx+1))
print(x)

cpos = np.zeros((Nt+1, Nx+1))
cneg = np.zeros((Nt+1, Nx+1))
rho = np.zeros((Nt+1, Nx+1))
rho_current = np.zeros(Nx+1)

for ix in np.arange(0, Nx):
    cpos[0,ix] = c_o
    cneg[0,ix] = c_o
    x[ix] = ix*dx
    #cpos[0,ix] = c_o*x[ix]*(Lx-x[ix])
    #cneg[0,ix] = c_o*x[ix]*(Lx-x[ix])
cmax = np.max(cpos)

#cpos[0,int(Nx/2-1)] = 1.0
# N_half = int(Nx/2)
# print(N_half)
#cpos[0,N_half] = 1.0
#cpos[0,int(.5/dx)] = 1.0

#plt.plot(x[:],cpos[0,:])

#matrix = np.zeros((Nx-1, Nx-1))
inverse_matrix = np.zeros((Nx+1, Nx+1))

#rho = np.zeros((Nx-1))
#rho[0] = 1.0

#phi = np.zeros((Nx))
phi = np.zeros((Nt+1, Nx+1))

phi_left = 5.0
phi_right = 0.0



rho[0,0] = rho_eo*dx**2 - phi_left 
rho[0,Nx] = rho_eo*dx**2 - phi_right

for i in np.arange(1, Nx-1):
    #rho[i] = 4.0*x[i]*(1.0-x[i])
    rho[0,i] = rho_eo*dx**2
#Calculation of Inverse 1D Poisson Matrix
#https://www.scirp.org/journal/paperinformation.aspx?paperid=50041    
for j in np.arange(0, gp_init):
    for i in np.arange(0, gp_init):
        if (i >= j):
            element = -(j+1)*(gp_init -i )/(gp_init+1)
        elif (i < j):
            element = -(i+1)*(gp_init -j )/(gp_init+1)
        inverse_matrix[i,j] = element

# print("I = ",np.shape(inverse_matrix))

rho_current = rho[0,0:Nx+1]
# print("rho_current = ",np.shape(rho_current) ) 
#print(rho_current)
phi[0,:] = inverse_matrix @ rho_current



for it in range(1, Nt+1) :
    print("it = ", it, " of ", Nt)
    for ix in np.arange(1, Nx) :
        #print(ix)
        if (ix == 1):
            cpos[it,ix] = cpos[it-1,ix] + (Dpos*(cpos[it-1,ix+1]-cpos[it-1,ix])/dx**2)*dt
            cpos[it,0] = cpos[it,1]
            cneg[it,ix] = cneg[it-1,ix] + (Dneg*(cneg[it-1,ix+1]-cneg[it-1,ix])/dx**2)*dt
            cneg[it,0] = cneg[it,1]
            #print("left")
        elif (ix == Nx-1) :
            cpos[it,ix] = cpos[it-1,ix] - (Dpos*(cpos[it-1,ix]-cpos[it-1,ix-1])/dx**2)*dt
            cpos[it,Nx] = cpos[it,Nx-1]
            cneg[it,ix] = cneg[it-1,ix] - (Dneg*(cneg[it-1,ix]-cneg[it-1,ix-1])/dx**2)*dt
            cneg[it,Nx] = cneg[it,Nx-1]
            #print("right")
        else :
            cpos[it,ix] = cpos[it-1,ix] + (Dpos*(cpos[it-1,ix+1]-2*cpos[it-1,ix]+cpos[it-1,ix-1])/dx**2)*dt
            cneg[it,ix] = cneg[it-1,ix] + (Dneg*(cneg[it-1,ix+1]-2*cneg[it-1,ix]+cneg[it-1,ix-1])/dx**2)*dt
            #print("mid")
        for i in np.arange(1, Nx-1):
            #rho[i] = 4.0*x[i]*(1.0-x[i])
            rho[it,i] = 1.602e-19*cpos[it,i]*dx**2
        rho[it,0] = 1.602e-19*cpos[it,0]*dx**2 - phi_left 
        rho[it,Nx] = 1.602e-19*cpos[it,Nx]*dx**2 - phi_right
        rho_current = rho[it,:]
        #print(rho_current)
        phi[it,:] = inverse_matrix @ rho_current
        

#fig1, ax = plt.subplots()

# def animate(i):
#     plt.xlim(0, 1.0)
#     time = i*dt
#     #ax.set_facecolor('#74ccf4')
#     plt.title("Concentration \n t = %5.2e s" %(time))
#     plt.plot(x[:],cpos[i,:], color='yellow')
#     #plt.plot(x[:],None)

# anim = animation.FuncAnimation(fig1, animate, init_func=None,frames=Nt, interval=0.1, blit=False, repeat=False)

# anim.save('image.gif', fps=120, writer="pillow")
#plt.ioff()
i = 0
# while (i <= Nt):
#     t_plot = i*dt
#     #fig, ax = plt.subplots()
#     plate = plt.figure(i, frameon=True, figsize=[8.0, 10.0],dpi = 300, constrained_layout= True)
#     #plt.tight_layout()
#     #ax.rcParams.update({'font.size': 18})
#     ax=plt.gca()   
#     #ax.rcParams['axes.facecolor'] = '#74ccf4'
#     #ax.set_facecolor('#74ccf4')
#     plt.xlim(0, Lx)
#     plt.ylim(0, 2.0*c_o )
#     plt.xticks(fontsize= 18)
#     plt.yticks(fontsize=18)
#     time = i*dt
#     #plt.tight_layout()
#     plt.title("Diffusion of NaCl in water \n t = %5.2f s/ %5.2f s" %(time, Ttotal), fontsize=18)
#     plt.plot(x,cpos[i,:], color='#ffc65a', label = "Na")
#     plt.plot(x,cneg[i,:], color = "#30583b", label = "Cl")
#     ax.yaxis.set_major_formatter(FormatStrFormatter('%5.2e'))
#     plt.xlabel("L [m]", fontsize=18)
#     # plt.ylabel("Concentration  [m^-3]", fontsize=12)
#     # r'\textbf{time (s)}'
#     #plt.ylabel(r'\textbf{Concentration}  $\mathbf{[m^{-3}]}$', fontsize=18)
#     #ax.set_ylabel(r'\textbf{Concentration}~  [m^{-3}]', fontsize=12)
#     plt.legend(loc="upper left", fontsize=18)
    

#     plate_name = "plate_%04d.png" % (i)
#     #plate.patch.set_alpha(0.0)
#     #plate.patch.set_alpha(1.0)

#     plt.savefig(plate_name, dpi='figure', format='png', bbox_inches=None, pad_inches=0.25, transparent=False)

#     i = i + 1
#     plt.ioff()
#     plt.close(plate)
#     #plt.show()


#cmd = "convert plate_*.png movie.gif"
cmd0 = "rm output.mp4"
os.system(cmd0)
cmd1 = "ffmpeg -framerate 10 -i plate_%04d.png output.mp4"
os.system(cmd1)
# cmd2 = "rm *.png"
# os.system(cmd2)

phi_plot = np.zeros((Nx + 1))
# phi_plot[0]= phi_left
# phi_plot[Nx] = phi_right

for i in np.arange(0, Nx+1):
    phi_plot[i] = phi[1,i]
    

# def phi_analytic(x) :
#     phi = 0.5*rho[0,x]*x**2 -(phi_left+rho_o/2.0)*x + phi_left
#     return phi

phi_plot_analytic = np.zeros((Nx+1))

C_analytic = (phi_right - phi_left - 0.5*rho_eo*Lx**2)/Lx
for i in np.arange(0, Nx+1):
    x[i] = i*dx
    phi_plot_analytic[i] = 0.5*rho_eo*x[i]**2 + C_analytic*x[i] + phi_left
#x[Nx] = Lx
print(cpos[0,0:Nx])
plt.ion()
plt.xlim(0, Lx)
plt.plot(x,phi_plot,label = "numerical")
plt.plot(x,phi_plot_analytic, label = "analytic")
plt.legend()
plt.show

#play sound when done
finished_sound = os.path.dirname(os.path.realpath(__file__)) + '/' + 'paddansq.wav'
play_done_sound = subprocess.run(["play","-v","0.25", finished_sound])





