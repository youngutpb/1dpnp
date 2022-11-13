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
Dneg = 2.03e-9# Cl- m^2/s

#Dpos = 1.0 #Na+ m^2/s
#Dneg = 1.0 # Cl- m^2/s


###############  Adjustable Parameters
Z_anion = 1.00
Z_cation = -1.0


plt.rcParams["font.family"] = "times"
plt.rcParams['text.usetex'] = True

Lx = 0.1
Nx = 100

Nt = 5
Nplot = 1

cn = 4.0E11
#cn = 0.0

phi_left = 0
phi_right = 1.5


#########################################################################

## Stability Criterion:  dt <= dx^2/(2*D)
## https://www.atmos.albany.edu/facstaff/brose/classes/ATM623_Spring2015/Notes/Lectures/Lecture16%20--%20Numerical%20methods%20for%20diffusion%20models.html


dx = Lx/Nx
#dt = 0.25
#dt_min = 2.0*np.trunc(0.5*(dx**2)/(2*max(Dpos,Dneg)))
dt_min = (dx**2)/(2*max(Dpos,Dneg))
print(dt_min)
#dt = dt_min
dt = dt_min
#dt = 2
Ttotal = Nt*dt
print(dt)
x = np.linspace(0,Lx, Nx+1)
#print(x)

cn_anion = np.zeros((Nt+1, Nx+1))
cn_cation = np.zeros((Nt+1, Nx+1))

matrix = np.zeros((Nx+1, Nx+1))
inverse = np.zeros((Nx+1, Nx+1))

rho = np.zeros((Nx+1))
rho_epo = np.zeros(Nx+1)

phi = np.zeros(( Nt+1,Nx+1))

#y=[]

#################### Calculation of Poisson Inverse Matrix ################33
for j in np.arange(0, Nx+1):
    for i in np.arange(0, Nx+1):
        if (i >= j):
            element = -(j+1)*(Nx+1-i)/(Nx+1+1)
        elif (i < j):
            element = -(i+1)*(Nx+1-j)/(Nx+1+1)
        inverse[i,j] = element


for ix in np.arange(0, Nx+1):
    #cn_anion[0,ix] = cn
    #cn_cation[0,ix] = cn
    # cn_anion[0,ix] = cn*x[ix]*(Lx-x[ix])
    # cn_cation[0,ix] = cn*x[ix]*(Lx-x[ix])
    #cn_anion[0,ix] = cn*x[ix]/Lx
    #cn_cation[0,ix] = cn*(1-x[ix]/Lx)
    #cn_anion[0,ix] = cn*np.exp(-x[ix]/(0.05*Lx))
    #cn_cation[0,ix] = cn*np.exp((x[ix]-Lx)/(0.05*Lx))
    n=10
    cn_anion[0,ix] = cn*((x[ix]-Lx)**n)/Lx**n
    cn_cation[0,ix] = cn*x[ix]**n/Lx**n
    
    # cn_cation[0,ix] = cn*((x[ix]-Lx)**n)/Lx**n
    # cn_anion[0,ix] = cn*x[ix]**n/Lx**n
    
    cmax = np.max(cn_cation)
    
    rho[ix] = (Z_anion*q_ele*cn_anion[0,ix]+Z_cation*q_ele*cn_cation[0,ix])
    if (ix ==0):
        rho_epo[0] = ((Z_anion*q_ele*cn_anion[0,ix]+Z_cation*q_ele*cn_cation[0,ix])/epo)*dx**2 - phi_left
    elif (ix==Nx):
        rho_epo[Nx] = ((Z_anion*q_ele*cn_anion[0,ix]+Z_cation*q_ele*cn_cation[0,ix]))/epo*dx**2 - phi_right  
    else:
        rho_epo[ix] = rho[ix]/epo*dx**2
       
phi[0,:] = inverse @ rho_epo




#################### Main Loop  ################################################33
for it in range(1, Nt) :
    #print("it = ", it, " of ", Nt+1)
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
            
       
#for i in np.arange(0, Nx+1):
    #rho[i] = 4.0*x[i]*(1.0-x[i])
    #rho[i] = rho_eo*dx**2    \
        rho[ix] = (Z_anion*q_ele*cn_anion[it,ix]+Z_cation*q_ele*cn_cation[it,ix])
        if (ix ==0):
            rho_epo[0] = ((Z_anion*q_ele*cn_anion[it,ix]+Z_cation*q_ele*cn_cation[it,ix])/epo)*dx**2 - phi_left
        elif (ix==Nx):
            rho_epo[Nx] = ((Z_anion*q_ele*cn_anion[it,ix]+Z_cation*q_ele*cn_cation[it,ix]))/epo*dx**2 - phi_right  
        else:
            rho_epo[ix] = rho[ix]/epo*dx**2
       
    phi[it,:] = inverse @ rho_epo

            
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

it = 0
iplot = 0


#ax.set_ylabel(r'\textbf{Concentration}~  [m^{-3}]', fontsize=12)

    
while (it <= Nt):
#     #ax.rcParams.update({'font.size': 18})

    if(it == 0):
      cn_plate = plt.figure(it, frameon=True, figsize=[3.0, 4.0],dpi = 300, constrained_layout= True)
      ax=plt.gca()   
      ax.set_facecolor('#74ccf4')
      # #ax.rcParams['axes.facecolor'] = '#74ccf4'
      # #ax.set_facecolor('#74ccf4')

      #plt.xticks(fontsize=18)
      #plt.yticks(fontsize=18)    
      plt.rcParams.update({'font.size': 12})
      plt.grid(visible=True, which="minor", axis='y', color='black', linestyle='-', linewidth=2)
      # #ax.yaxis.set_major_formatter(FormatStrFormatter('%5.2e'))
      plt.xlabel(r"\textbf{L~[m]}")
      plt.ylabel(r'\textbf{Concentration}  $\mathbf{[m^{-3}]}$')
      # #plt.tight_layout()
      # #ax.rcParams.update({'font.size': 18})
      # # ax=plt.gca()   
      # # #ax.rcParams['axes.facecolor'] = '#74ccf4'

      plt.xlim(0, Lx)
      plt.ylim(0, 5E13)

      time = it*dt/3600
      Ttotal_hr = Ttotal/3600

      plt.title("Diffusion of NaCl in water \n t = %5.2f hr/ %5.2f hr \n" %(time, Ttotal_hr), fontsize=12)
      
      plt.plot(x[:],cn_cation[it,:], color = "#30583b", label = "Cl")
      plt.plot(x,cn_anion[it,:], color='#ffc65a', label = "Na")

     
      plt.grid(visible=True, which="both", axis='both', color='black', linestyle='--', linewidth=1)         #plt.xlabel("L [m]", fontsize=18)
      # # plt.ylabel("Concentration  [m^-3]", fontsize=12)
      # # r'\textbf{time (s)}'
      plt.ylabel(r'\textbf{Concentration}  $\mathbf{[m^{-3}]}$', fontsize=12)

      # #ax.set_ylabel(r'\textbf{Concentration}~  [m^{-3}]', fontsize=12)
      plt.legend(loc="upper left", fontsize=12)
  
      # #plt.legend(loc="upper left", fontsize=18)
      cn_plate_name = "cn_plate_%04d.png" % (it/Nplot)
      # #plate.patch.set_alpha(0.0)
      # #plate.patch.set_alpha(1.0)

      plt.savefig(cn_plate_name, dpi='figure', format='png', bbox_inches=None, pad_inches=0.25, transparent=False)
#iplot = 0

    elif(iplot == Nplot):
        t_plot = i*dt
        # #fig, ax = plt.subplots()
        cn_plate = plt.figure(it, frameon=True, figsize=[3.0, 4.0],dpi = 300, constrained_layout= True)
        ax=plt.gca()   
        ax.set_facecolor('#74ccf4')
        # #ax.rcParams['axes.facecolor'] = '#74ccf4'
        # #ax.set_facecolor('#74ccf4')

        # plt.xticks(fontsize=18)
        # plt.yticks(fontsize=18)    
        plt.rcParams.update({'font.size': 12})
        plt.grid(visible=True, which="minor", axis='y', color='black', linestyle='-', linewidth=2)
        # #ax.yaxis.set_major_formatter(FormatStrFormatter('%5.2e'))
        plt.xlabel(r"\textbf{x~[m]}")
        plt.ylabel(r'\textbf{Concentration}  $\mathbf{[m^{-3}]}$')
        # #plt.tight_layout()
        # #ax.rcParams.update({'font.size': 18})
        # # ax=plt.gca()   
        # # #ax.rcParams['axes.facecolor'] = '#74ccf4'

        plt.xlim(0, Lx)
        plt.ylim(0, cmax)

        #time = it*dt
        time = it*dt/3600
        Ttotal_hr = Ttotal/3600
 
        plt.title("Diffusion of NaCl in water \n t = %5.2f hr/ %5.2f hr \n" %(time, Ttotal_hr), fontsize=12)
        
        plt.plot(x[:],cn_cation[it,:], color = "#30583b", label = "Cl")
        plt.plot(x,cn_anion[it,:], color='#ffc65a', label = "Na")

        #ax.yaxis.set_major_formatter(FormatStrFormatter('%5.2e'))
       
        plt.grid(visible=True, which="both", axis='both', color='black', linestyle='--', linewidth=1)         #plt.xlabel("L [m]", fontsize=18)
        # # plt.ylabel("Concentration  [m^-3]", fontsize=12)
        # # r'\textbf{time (s)}'
        plt.ylabel(r'\textbf{Concentration}  $\mathbf{[m^{-3}]}$', fontsize=12)

        # #ax.set_ylabel(r'\textbf{Concentration}~  [m^{-3}]', fontsize=12)
        plt.legend(loc="upper left", fontsize=12)
    
        # #plt.legend(loc="upper left", fontsize=18)
        cn_plate_name = "cn_plate_%04d.png" % (it/Nplot)
        # #plate.patch.set_alpha(0.0)
        # #plate.patch.set_alpha(1.0)

        plt.savefig(cn_plate_name, dpi='figure', format='png', bbox_inches=None, pad_inches=0.25, transparent=False)
        iplot = 0

    it = it + 1
    iplot = iplot +1
    plt.ioff()
    plt.close(cn_plate)
    


#print("################################################################")
#cmd = "convert plate_*.png movie.gif"
# cmd0 = "rm cn.mp4"
# os.system(cmd0)
# cmd1 = "ffmpeg -framerate 5 -i cn_plate_%04d.png cn.mp4"
# os.system(cmd1)
# cmd2 = "rm *.png"
# os.system(cmd2)



# for i in np.arange(1, Nx):
#     phi_plot[i] = phi[i-1]
    

########### Analytic Solution Computation for Check on Numerics ##############

# phi_analytic = np.zeros((Nx+1))

# for i in np.arange(0, Nx+1):
#     #x[i] = i*dx
#     ########### Solution for cn = constant
#     #C_analytic = phi_right/Lx - 0.5*rho[i]*Lx/epo - phi_left/Lx
#    # phi_analytic[i] = (0.5*rho[i]*x[i]**2)/epo + C_analytic*x[i] + phi_left
#    ########### Solution for cn_anion = cn*x/Lx, cn_cation = cn*(1-x/Lx)
#    #Wphi_analytic[i] = -((Z_cation-Z_anion)*cn*q_ele*x[i]**3-3*Lx*Z_cation*cn*q_ele*x[i]**2)/(6*Lx*epo) - (((2*Lx**2*Z_cation+Lx**2*Z_anion)*cn*q_ele-6*epo*phi_right+6*epo*phi_left)*x[i])/(6*Lx*epo) + phi_left
#    ############# Solution for cn_anion = cn_cation = cn(Lx-x)   ######################
#    phi_analytic[i] = -(Z_anion*cn*q_ele*x[i]**4+(2*Z_cation-2*Lx*Z_anion)*cn*q_ele*x[i]**3-6*Lx*Z_cation*cn*q_ele*x[i]**2)/(12*epo)- (((4*Lx**3*Z_cation+Lx**4*Z_anion)*cn*q_ele-12*epo*phi_right+12*epo*phi_left)*x[i])/(12*Lx*epo) - phi_left
    



# phi_plot = np.zeros((Nx + 1))
# phi_plot = phi[int(Nt/2),:]
# phi_plot[0]= phi_left
# phi_plot[Nx] = phi_right

it = 0
iplot = 0
while(it<=Nt):
    if(it ==0):
        #plt.ion()
        # t_plot = it*dt
        plt.rcParams.update({'font.size': 12})
        time = it*dt/3600
        Ttotal_hr = Ttotal/3600
        phi_plate = plt.figure(it, frameon=True, figsize=[3.0, 4.0],dpi = 300, constrained_layout= True)
        plt.xlim(0, Lx)
        plt.ylim(0, 5)
        plt.plot(x,phi[it,:], label = "numerical")
        plt.xlabel("x [m]")
        plt.ylabel("phi [volts]")
        plt.grid(visible=True, which="both", axis='both', color='black', linestyle='--', linewidth=1)         #plt.xlabel("L [m]", fontsize=18)
        phi_plate_name = "phi_plate_%04d.png" % (it)
        plt.title("Diffusion of NaCl in water \n t = %5.2f hr/ %5.2f hr \n $\phi_L = %3.1f~volts ~ \phi_R = %3.1f ~volts$ " %(time, Ttotal_hr, phi_left, phi_right), fontsize=12)
    #plate.patch.set_alpha(0.0)
    #plate.patch.set_alpha(1.0)

        plt.savefig(phi_plate_name, dpi='figure', format='png', bbox_inches=None, pad_inches=0.25, transparent=False)
    elif(iplot == Nplot):
        #plt.ion()
        # t_plot = it*dt
        plt.rcParams.update({'font.size': 12})
        #time = it*dt
        time = it*dt/3600
        Ttotal_hr = Ttotal/3600
        phi_plate = plt.figure(it, frameon=True, figsize=[3.0, 4.0],dpi = 300, constrained_layout= True)
        plt.xlim(0, Lx)
        #plt.ylim(0, 5)
        plt.plot(x,phi[it,:], label = "numerical")
        plt.xlabel("x [m]")
        plt.ylabel(r"$\phi~ [volts]$")
        plt.grid(visible=True, which="both", axis='both', color='black', linestyle='--', linewidth=1)         #plt.xlabel("L [m]", fontsize=18)
        phi_plate_name = "phi_plate_%04d.png" % (it/Nplot)
        plt.title("Diffusion of NaCl in water \n t = %5.2f hr/ %5.2f hr \n $\phi_L = 1.5~volts ~ \phi_R = 0.0~volts$ " %(time, Ttotal_hr), fontsize=12)
    #plate.patch.set_alpha(0.0)
    #plate.patch.set_alpha(1.0)

        plt.savefig(phi_plate_name, dpi='figure', format='png', bbox_inches=None, pad_inches=0.25, transparent=False)
        iplot = 0

    it = it + 1
    iplot = iplot +1
    plt.ioff()
    plt.close(phi_plate)
#plt.plot(x,phi_analytic, label = "analytic", linestyle = "--")
#plt.legend()
#plt.show

print("################################################################")

# cmd3 = "rm phi.mp4"
# os.system(cmd3)
# cmd4 = "ffmpeg -framerate 5 -i phi_plate_%04d.png phi.mp4"
# os.system(cmd4)
# cmd5 = "rm phi_*.png"
# os.system(cmd5)

#play sound when done
finished_sound = os.path.dirname(os.path.realpath(__file__)) + '/' + 'paddansq.wav'
play_done_sound = subprocess.run(["play","-v","0.25", finished_sound])





