#This is an attempt to make a program to simulate fluid flow using LBM
#This is for single-phase flow

####################################################### LIBRARIES ###########################################################################
import math
import numpy as np
from numpy import *
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import lbm_inputdataSCMP as ini             #import from theinput file

####################################################### SIMULATION ###########################################################################

for t in range(ini.T):

    #Zou and He velocity BCs on west side
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.m[i,j] == 2:
                ini.rho[i,j] = (ini.f[i,j,0] + ini.f[i,j,2] + ini.f[i,j,4] + 2.*(ini.f[i,j,3] + ini.f[i,j,7] + ini.f[i,j,6])) / (1 - ini.ux0)
                ru = ini.rho[i,j]*ini.ux0
                ini.f[i,j,1] = ini.f[i,j,3] + (2./3.)*ru
                ini.f[i,j,5] = ini.f[i,j,7] + (1./6.)*ru - (1./2.)*(ini.f[i,j,2] - ini.f[i,j,4])
                ini.f[i,j,8] = ini.f[i,j,6] + (1./6.)*ru - (1./2.)*(ini.f[i,j,4] - ini.f[i,j,2])
    
    #Zou and He velocity BCs on south side
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.m[i,j] == 3:
                ini.rho[i,j] = (ini.f[i,j,0] + ini.f[i,j,1] + ini.f[i,j,3] + 2.*(ini.f[i,j,4] + ini.f[i,j,7] + ini.f[i,j,8])) / (1 - ini.uy0)
                ru = ini.rho[i,j]*ini.uy0
                ini.f[i,j,2] = ini.f[i,j,4] + (2./3.)*ru
                ini.f[i,j,5] = ini.f[i,j,7] + (1./6.)*ru - (1./2.)*(ini.f[i,j,1] - ini.f[i,j,3])
                ini.f[i,j,6] = ini.f[i,j,8] + (1./6.)*ru - (1./2.)*(ini.f[i,j,3] - ini.f[i,j,1])
                

    # ... computing density for imaging
    ini.rho[:,:] = 0.
    for a in range(9):
        ini.f[:,:,a] = np.where(ini.f[:,:,a] > ini.f_tol, ini.f[:,:,a], ini.f_tol)
        ini.rho[:,:] += np.where(ini.m[i,j] == 0, ini.f[:,:,a], 0)
        
    #print ini.f

    if ini.f.min() < 0.:
        print "Operation terminated 1, negative probability found"
        print ini.f.min()
        break

    ###SCMP calculation

    #compute psi
    ini.psi[:,:] = where(ini.m[:,:] == 0, ini.psi0 * np.exp(-ini.rho0 / ini.rho[:,:]), 0)

    #compute pressure
    ini.P[:,:] = 0.
    ini.P[:,:] = ini.rho[:,:]/3 + ini.G*(ini.psi[:,:]**2)/6

    #compute interaction force
    for j in range(ini.sizeY_+2):
        j_n = (j-1) if j > 0 else (ini.sizeY_+1)
        j_p = (j+1) if j < (ini.sizeY_ + 1) else 0 
        for i in range(ini.sizeX_+ 2):
            i_n = (i-1) if i > 0 else (ini.sizeX_+1)
            i_p = (i+1) if i < (ini.sizeX_ + 1) else 0 

            ini.Fx[i,j] = 0.
            ini.Fy[i,j] = 0.

            if (ini.m[i,j] == 0):
                ini.Fx[i,j] += ini.w[1]*ini.e_[0,1]*ini.psi[i_p, j]
                ini.Fy[i,j] += ini.w[1]*ini.e_[1,1]*ini.psi[i_p, j]

                ini.Fx[i,j] += ini.w[1]*ini.e_[0,2]*ini.psi[i, j_p]
                ini.Fy[i,j] += ini.w[1]*ini.e_[1,2]*ini.psi[i, j_p]

                ini.Fx[i,j] += ini.w[1]*ini.e_[0,3]*ini.psi[i_n, j]
                ini.Fy[i,j] += ini.w[1]*ini.e_[1,3]*ini.psi[i_n, j]

                ini.Fx[i,j] += ini.w[1]*ini.e_[0,4]*ini.psi[i, j_n]
                ini.Fy[i,j] += ini.w[1]*ini.e_[1,4]*ini.psi[i, j_n]

                ini.Fx[i,j] += ini.w[2]*ini.e_[0,1]*ini.psi[i_p, j_p]
                ini.Fy[i,j] += ini.w[2]*ini.e_[1,1]*ini.psi[i_p, j_p]

                ini.Fx[i,j] += ini.w[2]*ini.e_[0,2]*ini.psi[i_n, j_p]
                ini.Fy[i,j] += ini.w[2]*ini.e_[1,2]*ini.psi[i_n, j_p]

                ini.Fx[i,j] += ini.w[2]*ini.e_[0,3]*ini.psi[i_n, j_n]
                ini.Fy[i,j] += ini.w[2]*ini.e_[1,3]*ini.psi[i_n, j_n]

                ini.Fx[i,j] += ini.w[2]*ini.e_[0,4]*ini.psi[i_p, j_n]
                ini.Fy[i,j] += ini.w[2]*ini.e_[1,4]*ini.psi[i_p, j_n]

                ini.Fx[i,j] = -ini.G*ini.psi[i,j] * ini.Fx[i,j]
                ini.Fy[i,j] = -ini.G*ini.psi[i,j] * ini.Fy[i,j]

    #Streaming step
    for j in range(ini.sizeY_+2):
        j_n = (j-1) if j > 0 else (ini.sizeY_+1)
        j_p = (j+1) if j < (ini.sizeY_ + 1) else 0 
        for i in range(ini.sizeX_+ 2):
            i_n = (i-1) if i > 0 else (ini.sizeX_+1)
            i_p = (i+1) if i < (ini.sizeX_ + 1) else 0 
            if (ini.m[i,j] <> 1):
                #For streaming part, if adjacent grid is a m, then the density distribution will propagate, 
                #Else, the density is bounced back to the same grid, but with different direction
                
                ini.ftemp[i,j,0] = ini.f[i,j,0]
                
                if (ini.m[i_p,j] <> 1):    ini.ftemp[i_p,j,1] = ini.f[i,j,1]
                else:                           ini.ftemp[i,j,3] = ini.f[i,j,1]

                if (ini.m[i,j_p] <> 1):    ini.ftemp[i][j_p][2] = ini.f[i][j][2]
                else:                           ini.ftemp[i][j][4] = ini.f[i][j][2]

                if (ini.m[i_n,j] <> 1):    ini.ftemp[i_n,j,3] = ini.f[i,j,3]
                else:                           ini.ftemp[i,j,1] = ini.f[i,j,3]

                if (ini.m[i,j_n] <> 1):    ini.ftemp[i,j_n,4] = ini.f[i,j,4]
                else:                           ini.ftemp[i,j,2] = ini.f[i,j,4]

                if (ini.m[i_p,j_p] <> 1):  ini.ftemp[i_p,j_p,5] = ini.f[i,j,5]
                else:                           ini.ftemp[i,j,7] = ini.f[i,j,5]

                if (ini.m[i_n,j_p] <> 1):  ini.ftemp[i_n,j_p,6] = ini.f[i,j,6]
                else:                           ini.ftemp[i,j,8] = ini.f[i,j,6]

                if (ini.m[i_n,j_n] <> 1):  ini.ftemp[i_n,j_n,7] = ini.f[i,j,7]
                else:                           ini.ftemp[i,j,5] = ini.f[i,j,7]

                if (ini.m[i_p,j_n] <> 1):  ini.ftemp[i_p,j_n,8] = ini.f[i,j,8]
                else:                           ini.ftemp[i,j,6] = ini.f[i,j,8]
    
    # ... and then computing macroscopic density and velocity for each lattice point, after shifting
    ini.rho[:,:] = 0.
    ini.ux[:,:] = 0.
    ini.uy[:,:] = 0.
    for a in range(9):
        ini.rho[:,:] += ini.ftemp[:,:,a]    
        ini.ux[:,:] += ini.e_[0,a]*ini.ftemp[:,:,a]
        ini.uy[:,:] += ini.e_[1,a]*ini.ftemp[:,:,a]
    ini.ux[:,:] = np.where(ini.rho[:,:] > ini.f_tol, ini.ux[:,:]/ini.rho[:,:], ini.ux[:,:]/ini.rho[:,:])
    ini.uy[:,:] = np.where(ini.rho[:,:] > ini.f_tol, ini.uy[:,:]/ini.rho[:,:], ini.uy[:,:]/ini.rho[:,:])
    ini.u[:,:] = sqrt((ini.ux[:,:]**2 + ini.uy[:,:]**2)/2)

    #Computing equilibrium distribution function
    for j in range(ini.sizeY_+2):
        for i in range(ini.sizeX_+ 2):
            if ini.m[i,j] == 0:
                fct1 = ini.w[0]*ini.rho[i,j]
                fct2 = ini.w[1]*ini.rho[i,j]
                fct3 = ini.w[2]*ini.rho[i,j]

                ini.uxeq[i,j] = ini.ux[i,j]  + ini.tau0*ini.Fx[i,j]/ini.rho[i,j]                                     #uxeq will incorporate external forces, if any
                ini.uyeq[i,j] = ini.uy[i,j]  + ini.tau0*ini.Fy[i,j]/ini.rho[i,j]

                uxsq = ini.uxeq[i,j]*ini.uxeq[i,j]
                uysq = ini.uyeq[i,j]*ini.uyeq[i,j]

                uxuy5 = ini.uxeq[i,j] + ini.uyeq[i,j]
                uxuy6 = -ini.uxeq[i,j] + ini.uyeq[i,j]
                uxuy7 = -ini.uxeq[i,j] - ini.uyeq[i,j]
                uxuy8 = ini.uxeq[i,j] - ini.uyeq[i,j]

                usq = uxsq + uysq

                ini.feq[i,j,0] = fct1*(1.                                                   - ini.c_eq[2]*usq)
                ini.feq[i,j,1] = fct2*(1. + ini.c_eq[0]*ini.uxeq[i,j] + ini.c_eq[1]*uxsq    - ini.c_eq[2]*usq)
                ini.feq[i,j,2] = fct2*(1. + ini.c_eq[0]*ini.uyeq[i,j] + ini.c_eq[1]*uysq    - ini.c_eq[2]*usq)
                ini.feq[i,j,3] = fct2*(1. - ini.c_eq[0]*ini.uxeq[i,j] + ini.c_eq[1]*uxsq    - ini.c_eq[2]*usq)
                ini.feq[i,j,4] = fct2*(1. - ini.c_eq[0]*ini.uyeq[i,j] + ini.c_eq[1]*uysq    - ini.c_eq[2]*usq)
                ini.feq[i,j,5] = fct3*(1. + ini.c_eq[0]*uxuy5 + ini.c_eq[1]*uxuy5*uxuy5  - ini.c_eq[2]*usq)
                ini.feq[i,j,6] = fct3*(1. + ini.c_eq[0]*uxuy6 + ini.c_eq[1]*uxuy6*uxuy6  - ini.c_eq[2]*usq)
                ini.feq[i,j,7] = fct3*(1. + ini.c_eq[0]*uxuy7 + ini.c_eq[1]*uxuy7*uxuy7  - ini.c_eq[2]*usq)
                ini.feq[i,j,8] = fct3*(1. + ini.c_eq[0]*uxuy8 + ini.c_eq[1]*uxuy8*uxuy8  - ini.c_eq[2]*usq)

    #Collision step

    # for j in range(ini.sizeY_+2):
    #     for i in range(ini.sizeX_+ 2):
    #         if ini.m[i,j] == 0:
    #             for a in range(9):
    #                 ini.tau[i,j,a] = maximum(ini.tau0, (ini.ftemp[i,j,a] - ini.feq[i,j,a])/(ini.ftemp[i,j,a] + ini.f_tol))
    ini.tau[:,:,:] = maximum(ini.tau0, (1 - (ini.feq[:,:,:]/ini.ftemp[:,:,:])))

    print ini.tau
    for a in range(9):
       ini.f[:,:,a] = np.where(ini.m[:,:] == 0, ini.ftemp[:,:,a] - ((ini.ftemp[:,:,a] - ini.feq[:,:,a]) / ini.tau.max()), ini.f[:,:,a])

    #print ini.f

    if ini.f.min() < -0.000001:
        print "Operation terminated 2, negative probability found"
        print ini.f.min()
        break

    #Plotting 
    print "Time = ", t
    print "Mass = ", sum(ini.rho)
    print "Velocity x dir = ", sum(ini.ux)
    print "Velocity y dir = ", sum(ini.uy)
    print ini.f.min()

    if mod(t,1) == 0:
        varm = ini.P.transpose()        #Change the variable to the one that will be plotted: rho, ux, or uy
        plt.figure(1)
        #cmap = sns.cm.rocket_r
        #ax = sns.heatmap(varm, annot=False, vmin=0., vmax=0.6, cmap='RdYlBu_r')           #For ux
        ax = sns.heatmap(varm, annot=False, vmin=0, vmax=30, cmap='RdYlBu_r')           #For rho
        ax.invert_yaxis()
        plt.pause(0.001)
        plt.clf()


    #Plot cross-section velocity
    # plt.figure(2)
    # plt.plot(ini.ux[-2,:],label='Outlet')
    # plt.plot(ini.ux[(ini.sizeX_+2)//2,:],label='Middle')
    # plt.plot(ini.ux[2,:],label='Inlet')
    # plt.legend(loc='upper right')
    # plt.show
    # plt.pause(0.001)
    # plt.clf()

####################################################### OUTPUT ###########################################################################

        
