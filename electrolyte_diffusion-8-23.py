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
Lx = 1.0e-2
dx = Lx/Nx


# dx = 1.0e-5
# Lx = 0.2
# Nx = int(Lx/dx)

Nt = 4

# https://www.aqion.de/site/diffusion-coefficients
# https://www.physiologyweb.com/calculators/diffusion_time_calculator.html
Dpos = 1.33e-9 #Na+ m^2/s
Dneg = 2.03e-9 # Cl- m^2/s


## Stability Criterion:  dt <= dx^2/(2*D)
## https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture16%20--%20Numerical%20methods%20for%20diffusion%20models.html


#dt = 0.25
dt_min = dx**2/(2*max(Dpos,Dneg))
print(dt_min)
dt = dt_min
Ttotal = Nt*dt

x = np.linspace(0,Lx, Nx+1)
print(x)

cpos = np.zeros((Nt+1, Nx+1))
cneg = np.zeros((Nt+1, Nx+1))

for ix in np.arange(0, Nx):
    # cpos[0,ix] = 1.0
    # cneg[0,ix] = 1.0
    cpos[0,ix] = 4.0*x[ix]*(Lx-x[ix])
    cneg[0,ix] = 4.0*x[ix]*(Lx-x[ix])
    cmax = np.max(cneg)

#cpos[0,int(Nx/2-1)] = 1.0
# N_half = int(Nx/2)
# print(N_half)
#cpos[0,N_half] = 1.0
#cpos[0,int(.5/dx)] = 1.0

#plt.plot(x[:],cpos[0,:])



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
    fig, ax = plt.subplots()
    plate = plt.figure(i, frameon=True, figsize=[6.0, 6.0])
    #plt.tight_layout()
    #ax.rcParams.update({'font.size': 18})
    ax=plt.gca()   
    ax.rcParams['axes.facecolor'] = '#74ccf4'
    ax.plt.xlim(0, Lx)
    ax.plt.ylim(0, cmax)
    time = i*dt
    #plt.tight_layout()
    ax.plt.title("Concentration \n t = %5.2f s/ %5.2f s" %(time, Ttotal))
    ax.plt.plot(x[:],cpos[i,:], color='#ffc65a', label = "Na")
    ax.plt.plot(x[:],cneg[i,:], color = "#30583b", label = "Cl")
    #plt.legend()
    ax.plt.xlabel(" L [m]", fontsize=12)
    ax.plt.ylabel(" Concentration  [m^-3]", fontsize=12)
    ax.plt.legend(loc="upper left")
    

    plate_name = "plate_%04d.png" % (i)
    #plate.patch.set_alpha(0.0)
    #plate.patch.set_alpha(1.0)

    ax.plt.savefig(plate_name, dpi='figure', format='png', bbox_inches=None, pad_inches=0.1, transparent=False)

    i = i + 1
    ax.plt.ioff()
    ax.plt.close(plate)
    #plt.show()


cmd = "convert plate_*.png movie.gif"
os.system(cmd)
cmd2 = "rm *.png"
os.system(cmd2)

#play sound when done
finished_sound = os.path.dirname(os.path.realpath(__file__)) + '/' + 'paddansq.wav'
play_done_sound = subprocess.run(["play","-v","0.25", finished_sound])





