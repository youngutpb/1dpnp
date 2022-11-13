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

#############  Physical Constants ################################
q_ele = 1.602E-19
epo = 8.854E-12

# https://www.aqion.de/site/diffusion-coefficients
# https://www.physiologyweb.com/calculators/diffusion_time_calculator.html
Dpos = 1.33e-9 #Na+ m^2/s
Dneg = 2.03e-9 # Cl- m^2/s


###############  Adjustable Parameters
Z_anion = 1.00
Z_cation = -1.0


plt.rcParams["font.family"] = "times"
plt.rcParams['text.usetex'] = True

Lx = 1.0e-1
Nx = 10

Nt = 200

cn = 4.0E13
#cn = 0.0

phi_left = 5.0
phi_right = 0.0


#########################################################################

## Stability Criterion:  dt <= dx^2/(2*D)
## https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture16%20--%20Numerical%20methods%20for%20diffusion%20models.html


dx = Lx/Nx
#dt = 0.25
dt_min = dx**2/(2*max(Dpos,Dneg))
#print(dt_min)
dt = dt_min
Ttotal = Nt*dt

x = np.linspace(0,Lx, Nx+1)
#print(x)

cn_anion = np.zeros((Nt+1, Nx+1))
cn_cation = np.zeros((Nt+1, Nx+1))

for ix in np.arange(0, Nx):
    #cn_anion[0,ix] = cn
    #cn_cation[0,ix] = cn
    #cn_anion[0,ix] = cn*x[ix]*(Lx-x[ix])/Lx
    #cn_cation[0,ix] = cn*x[ix]*(Lx-x[ix])/Lx
    
    cn_anion[0,ix] = cn*((x[ix]-Lx)**2)/Lx**2
    cn_cation[0,ix] = cn*x[ix]**2/Lx**2
    
    #cn_anion[0,ix] = cn*x[ix]/Lx
    #cn_cation[0,ix] = cn*(1-x[ix]/Lx)
    cmax = np.max(cn_cation)
    
#rho_eo = cn*1.602e-19/8.854e-12
#rho_eo = cn*q_ele/epo

#cn_anion[0,int(Nx/2-1)] = 1.0
# N_half = int(Nx/2)
# print(N_half)
#cn_anion[0,N_half] = 1.0
#cn_anion[0,int(.5/dx)] = 1.0

#plt.plot(x[:],cn_anion[0,:])

matrix = np.zeros((Nx+1, Nx+1))
inverse = np.zeros((Nx+1, Nx+1))

rho = np.zeros((Nx+1))
rho_epo = np.zeros(Nx+1)

phi = np.zeros((Nx+1))


for i in np.arange(0, Nx+1):
    #rho[i] = 4.0*x[i]*(1.0-x[i])
    #rho[i] = rho_eo*dx**2    \
    rho[i] = (Z_anion*q_ele*cn_anion[0,i]+Z_cation*q_ele*cn_cation[0,i])
    if (i ==0):
        rho_epo[0] = ((Z_anion*q_ele*cn_anion[0,i]+Z_cation*q_ele*cn_cation[0,i])/epo)*dx**2 - phi_left
    elif (i==Nx):
       rho_epo[Nx] = ((Z_anion*q_ele*cn_anion[0,i]+Z_cation*q_ele*cn_cation[0,i]))/epo*dx**2 - phi_right  
    else:
        rho_epo[i] = rho[i]/epo*dx**2
       

print("rho =", rho)
    
for j in np.arange(0, Nx+1):
    for i in np.arange(0, Nx+1):
        if (i >= j):
            element = -(j+1)*(Nx+1-i)/(Nx+1+1)
        elif (i < j):
            element = -(i+1)*(Nx+1-j)/(Nx+1+1)
        inverse[i,j] = element
    
        
phi = inverse @ rho_epo



for it in range(1, Nt) :
    print("it = ", it, " of ", Nt)
    for ix in np.arange(1, Nx) :
        #print(ix)
        if (ix == 1):
            cn_anion[it,ix] = cn_anion[it-1,ix] + (Dpos*(cn_anion[it-1,ix+1]-cn_anion[it-1,ix])/dx**2)*dt
            cn_anion[it,0] = cn_anion[it,1]
            cn_cation[it,ix] = cn_cation[it-1,ix] + (Dneg*(cn_cation[it-1,ix+1]-cn_cation[it-1,ix])/dx**2)*dt
            cn_cation[it,0] = cn_cation[it,1]
            #print("left")
        elif (ix == Nx-1) :
            cn_anion[it,ix] = cn_anion[it-1,ix] - (Dpos*(cn_anion[it-1,ix]-cn_anion[it-1,ix-1])/dx**2)*dt
            cn_anion[it,Nx] = cn_anion[it,Nx-1]
            cn_cation[it,ix] = cn_cation[it-1,ix] - (Dneg*(cn_cation[it-1,ix]-cn_cation[it-1,ix-1])/dx**2)*dt
            cn_cation[it,Nx] = cn_cation[it,Nx-1]
            #print("right")
        else :
            cn_anion[it,ix] = cn_anion[it-1,ix] + (Dpos*(cn_anion[it-1,ix+1]-2*cn_anion[it-1,ix]+cn_anion[it-1,ix-1])/dx**2)*dt
            cn_cation[it,ix] = cn_cation[it-1,ix] + (Dneg*(cn_cation[it-1,ix+1]-2*cn_cation[it-1,ix]+cn_cation[it-1,ix-1])/dx**2)*dt
            #print("mid")
#     plt.plot(x[:],cn_anion[it,:])
# plt.show()

#fig1, ax = plt.subplots()

# def animate(i):
#     plt.xlim(0, 1.0)
#     time = i*dt
#     #ax.set_facecolor('#74ccf4')
#     plt.title("Concentration \n t = %5.2e s" %(time))
#     plt.plot(x[:],cn_anion[i,:], color='yellow')
#     #plt.plot(x[:],None)

# anim = animation.FuncAnimation(fig1, animate, init_func=None,frames=Nt, interval=0.1, blit=False, repeat=False)

# anim.save('image.gif', fps=120, writer="pillow")

i = 0
while (i <= Nt):
    t_plot = i*dt
    #fig, ax = plt.subplots()
    plate = plt.figure(i, frameon=True, figsize=[8.0, 10.0],dpi = 300, constrained_layout= True)
    #plt.tight_layout()
    #ax.rcParams.update({'font.size': 18})
    ax=plt.gca()   
    #ax.rcParams['axes.facecolor'] = '#74ccf4'
    ax.set_facecolor('#74ccf4')
    plt.xlim(0, Lx)
    plt.ylim(0, cmax)
    plt.xticks(fontsize= 18)
    plt.yticks(fontsize=18)
    time = i*dt
    #plt.tight_layout()
    plt.title("Diffusion of NaCl in water \n t = %5.2f s/ %5.2f s" %(time, Ttotal), fontsize=18)
    plt.plot(x[:],cn_anion[i,:], color='#ffc65a', label = "Na")
    plt.plot(x[:],cn_cation[i,:], color = "#30583b", label = "Cl")
    ax.yaxis.set_major_formatter(FormatStrFormatter('%5.2e'))
    plt.xlabel("L [m]", fontsize=18)
    # plt.ylabel("Concentration  [m^-3]", fontsize=12)
    # r'\textbf{time (s)}'
    plt.ylabel(r'\textbf{Concentration}  $\mathbf{[m^{-3}]}$', fontsize=18)
    #ax.set_ylabel(r'\textbf{Concentration}~  [m^{-3}]', fontsize=12)
    plt.legend(loc="upper left", fontsize=18)
    

    plate_name = "plate_%04d.png" % (i)
    #plate.patch.set_alpha(0.0)
    #plate.patch.set_alpha(1.0)

    plt.savefig(plate_name, dpi='figure', format='png', bbox_inches=None, pad_inches=0.25, transparent=False)

    i = i + 1
    plt.ioff()
    plt.close(plate)
    #plt.show()


#cmd = "convert plate_*.png movie.gif"
cmd0 = "rm output.mp4"
os.system(cmd0)
cmd1 = "ffmpeg -framerate 10 -i plate_%04d.png output.mp4"
os.system(cmd1)
cmd2 = "rm *.png"
os.system(cmd2)

# phi_plot = np.zeros((Nx + 1))
# phi_plot[0]= phi_left
# phi_plot[Nx] = phi_right

# for i in np.arange(1, Nx):
#     phi_plot[i] = phi[i-1]
    

phi_analytic = np.zeros((Nx+1))

for i in np.arange(0, Nx+1):
    #x[i] = i*dx
    ########### Solution for cn = constant
    #C_analytic = phi_right/Lx - 0.5*rho[i]*Lx/epo - phi_left/Lx
   # phi_analytic[i] = (0.5*rho[i]*x[i]**2)/epo + C_analytic*x[i] + phi_left
   ########### Solution for cn_anion = cn*x/Lx, cn_cation = cn*(1-x/Lx)
   #Wphi_analytic[i] = -((Z_cation-Z_anion)*cn*q_ele*x[i]**3-3*Lx*Z_cation*cn*q_ele*x[i]**2)/(6*Lx*epo) - (((2*Lx**2*Z_cation+Lx**2*Z_anion)*cn*q_ele-6*epo*phi_right+6*epo*phi_left)*x[i])/(6*Lx*epo) + phi_left
   ############# Solution for cn_anion = cn_cation = cn*x*(Lx-x)/Lx   ######################
   #phi_analytic[i] = -(Z_anion*cn*q_ele*x[i]**4+(2*Z_cation-2*Lx*Z_anion)*cn*q_ele*x[i]**3-6*Lx*Z_cation*cn*q_ele*x[i]**2)/(12*epo)- (((4*Lx**3*Z_cation+Lx**4*Z_anion)*cn*q_ele-12*epo*phi_right+12*epo*phi_left)*x[i])/(12*Lx*epo) - phi_left
   phi_analytic[i] = ((Z_cation+Z_anion)*cn*q_ele*x[i]**2)/(2*epo) - (((Lx**2*Z_cation+Lx**2*Z_anion)*cn*q_ele-2*epo*phi_right+2*epo*phi_left)*x[i])/(2*Lx*epo) + phi_left
 


plt.ion()
plt.xlim(0, Lx)
plt.plot(x,phi, label = "numerical")
plt.plot(x,phi_analytic, label = "analytic", linestyle = "--")
plt.legend()
plt.show

#play sound when done
finished_sound = os.path.dirname(os.path.realpath(__file__)) + '/' + 'paddansq.wav'
play_done_sound = subprocess.run(["play","-v","0.25", finished_sound])





