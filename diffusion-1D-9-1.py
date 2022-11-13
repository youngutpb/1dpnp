#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 16 23:55:21 2022

@author: douglas
"""
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.animation as animation
import os


Nx = 100
Lx = 0.2
dx = Lx/Nx
#dx = 0.1e-2
#Nx = int(Lx/dx)
#D = 0.01
D = 1.33e-5 

## Stability Criterion:  dt <= dx^2/(2*D)
## https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture16%20--%20Numerical%20methods%20for%20diffusion%20models.html

#Nt = 40
#dt = 0.25
dt_min = dx**2/(2*D)
print(dt_min)
dt = min(dt_min, 0.1)
#Ttotal = Nt*dt
Ttotal = 120
Nt = int(Ttotal/dt)
print("Nt = ", Nt)
x = np.linspace(0,Lx, Nx+1)
print(x)

c = np.zeros((Nt+1, Nx+1))

for ix in np.arange(0, Nx):
    c[0,ix] = 4.0*x[ix]*(Lx-x[ix]) 
    cmax = np.max(c)

#c[0,int(Nx/2-1)] = 1.0
# N_half = int(Nx/2)
# print(N_half)
#c[0,N_half] = 1.0
#c[0,int(.5/dx)] = 1.0

#plt.plot(x[:],c[0,:])



for it in range(1, Nt) :
    print("it = ", it)
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
#     plt.plot(x[:],c[it,:])
# plt.show()

#fig1, ax = plt.subplots()

# def animate(i):
#     plt.xlim(0, 1.0)
#     time = i*dt
#     #ax.set_facecolor('#74ccf4')
#     plt.title("Concentration \n t = %5.2e s" %(time))
#     plt.plot(x[:],c[i,:], color='yellow')
#     #plt.plot(x[:],None)

# anim = animation.FuncAnimation(fig1, animate, init_func=None,frames=Nt, interval=0.1, blit=False, repeat=False)

# anim.save('image.gif', fps=120, writer="pillow")

i = 1
while (i <= Nt):
    t_plot = i*dt
    plate = plt.figure(i, frameon=True, figsize=[6.0, 6.0])
    plt.xlim(0, Lx)
    plt.ylim(0, cmax)
    plt.tight_layout()
    time = i*dt
    plt.rcParams.update({'font.size': 18})
    axes=plt.gca()   
    plt.rcParams['axes.facecolor'] = '#74ccf4'
    plt.title("Concentration \n t = %5.2f s/%5.2f s" %(time, Ttotal))
    plt.plot(x[:],c[i,:], color='yellow')
    
    # plt.rcParams.update({'font.size': 14})
    # axes = plt.gca()
    # axes.set_aspect(1)
    # plt.rcParams['axes.facecolor'] = '#74ccf4'
    # #plt.rcParams['axes.facecolor'] = '#FFFFF0'
    # #plt.plot(x_Na,y_Na, color = ion_color_Na, linewidth=3, label = "Na")
    # #plt.plot(x_Cl,y_Cl, color = ion_color_Cl, linewidth=3, label = "Cl")
    # plt.scatter(x_Na[i-1], y_Na[i-1], c=ion_color_Na,
    #             marker="o", label="Na")
    # plt.scatter(x_Cl[i-1], y_Cl[i-1], c=ion_color_Cl,
    #             marker="o", label="Cl")
    # x_vec1, y_vec1 = np.meshgrid(
    #     np.linspace(-L/2., L/2, 10), np.linspace(0, H, 5))
    # x_vec2, y_vec2 = np.meshgrid(
    #     np.linspace(-L/2., L/2, 10), np.linspace(0.025, H-0.025, 4))
    # plt.scatter(x_vec1, y_vec1, 512, c='crimson', alpha=0.5,
    #             marker='$\u2190$', label="Electric Field", )
    # plt.scatter(x_vec2, y_vec2, 128, c='purple',
    #             marker='$\u29BF$', label="Magnetic Field", )
    # plt.legend('', frameon=False)
    # plt.legend(facecolor='white', framealpha=1)
    # plt.xlim(-L/2., L/2.)
    # plt.ylim(0.0, H)
    # plt.xlabel('x [m]', loc='center')
    # plt.ylabel('y [m]', loc='center')
    # plt.title('Position \n  Time =   %5.2e  s' %
    #           (t_plot), backgroundcolor='silver')
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





