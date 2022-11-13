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


Nx = 40
dx = 1.0/Nx
Nt = 200

gp_int = Nx - 1
rho_o = 1.0/8.854e-12

Dpos = 0.01
Dneg = 0.015


## Stability Criterion:  dt <= dx^2/(2*D)
## https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture16%20--%20Numerical%20methods%20for%20diffusion%20models.html


#dt = 0.25
dt_min = dx**2/(2*max(Dpos,Dneg))
print(dt_min)
dt = dt_min
Ttotal = Nt*dt

x = np.linspace(0,1.0, Nx+1)
#print(x)

cpos = np.zeros((Nt+1, Nx+1))
cneg = np.zeros((Nt+1, Nx+1))

matrix = np.zeros((Nx-1, Nx-1))
inverse_matrix = np.zeros((Nx-1, Nx-1))

rho = np.zeros((Nx-1))
#rho[0] = 1.0

phi = np.zeros((Nx))
Ex = np.zeros((Nx))

phi_left = 5.0
phi_right = 0.0



rho[0] = rho_o*dx**2 - phi_left
rho[Nx-2] = rho_o*dx**2 - phi_right 

for i in np.arange(1, Nx -2):
    #rho[i] = 4.0*x[i]*(1.0-x[i])
    rho[i] = rho_o*dx**2


for ix in np.arange(0, Nx):
    # cpos[0,ix] = 1.0
    # cneg[0,ix] = 1.0
    cpos[0,ix] = 4.0*x[ix]*(1.0-x[ix]) 
    cneg[0,ix] = 4.0*x[ix]*(1.0-x[ix]) 

#cpos[0,int(Nx/2-1)] = 1.0
# N_half = int(Nx/2)
# print(N_half)
#cpos[0,N_half] = 1.0
#cpos[0,int(.5/dx)] = 1.0

#plt.plot(x[:],cpos[0,:])
for j in np.arange(0,gp_int):
    for i in np.arange(0,gp_int):
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


for it in range(1, Nt) :
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
    phi = inverse_matrix @ rho
    
    
    
            #print("mid")
#     plt.plot(x[:],cpos[it,:])
# plt.show()

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

i = 0
while (i <= Nt):
    t_plot = i*dt
    plate = plt.figure(i, frameon=True, figsize=[6.0, 6.0])
    plt.tight_layout()
    plt.rcParams.update({'font.size': 18})
    axes=plt.gca()   
    plt.rcParams['axes.facecolor'] = '#74ccf4'
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.25)
    time = i*dt
    #plt.tight_layout()
    plt.title("Concentration \n t = %5.2e s" %(time))
    plt.plot(x[:],cpos[i,:], color='#ffc65a')
    plt.plot(x[:],cneg[i,:], color = "#30583b")
    

    plate_name = "plate_%04d.png" % (i)
    #plate.patch.set_alpha(0.0)
    #plate.patch.set_alpha(1.0)

    plt.savefig(plate_name, dpi='figure', format='png', bbox_inches=None, pad_inches=0.1, transparent=False)

    i = i + 1
    plt.ioff()
    plt.close(plate)
    #plt.show()


cmd = "convert plate_*.png movie.gif"
os.system(cmd)
cmd2 = "rm *.png"
os.system(cmd2)

#play sound when done
finished_sound = os.path.dirname(os.path.realpath(__file__)) + '/' + 'paddansq.wav'
play_done_sound = subprocess.run(["play","-v","0.25", finished_sound])





