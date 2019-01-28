#This is an attempt to make a program to simulate fluid flow using LBM
#This is for single-phase flow

####################################################### LIBRARIES ###########################################################################
import math
import numpy as np
from numpy import *
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import lbm_inputdata_klinkenberg as ini             #import from theinput file
import os

####################################################### SIMULATION ###########################################################################

# #... computing density for imaging
# ini.rho[:,:] = 0.
# for a in range(9):
#     ini.f[:,:,a] = np.where(abs(ini.f[:,:,a]) > ini.f_tol, ini.f[:,:,a], ini.f_tol)
#     ini.rho[:,:] += np.where(ini.m[:,:] <> 1, ini.f[:,:,a], 0)
# ini.rho[:,:] = ini.rho[:,:]*ini.rho0

for t in range(ini.T):

    #For smoothing, test friday 19/10/2018
    #if t == 0:
        #ini.tau0 = 2.*ini.tau0

    #if t == 20:
        #ini.tau0 = ini.tau0/2.

    # #west part
    # ini.rho[1,:] = (ini.f[1,:,0] + ini.f[1,:,2] + ini.f[1,:,4] + 2.*(ini.f[1,:,3] + ini.f[1,:,7] + ini.f[1,:,6])) / (1 - ini.ux0*(1 + 1e-4*sin(2.*arange(ini.sizeY_+2)*pi/ini.sizeY_)))
    # ru = ini.rho[1,:]*ini.ux0*(1 + 1e-4*sin(2*arange(ini.sizeY_+2)*pi/ini.sizeY_ + 2))
    # ini.f[1,:,1] = abs(ini.f[1,:,3] + (2./3.)*ru)
    # ini.f[1,:,5] = abs(ini.f[1,:,7] + (1./6.)*ru - (1./2.)*(ini.f[1,:,2] - ini.f[1,:,4]))
    # ini.f[1,:,8] = abs(ini.f[1,:,6] + (1./6.)*ru - (1./2.)*(ini.f[1,:,4] - ini.f[1,:,2]))

    #Pressure BC on the left
    #OK
    ini.ru[:] = ini.rho0*(ini.rho_left - ((ini.ftemp[1,:,0] + ini.ftemp[1,:,2] + ini.ftemp[1,:,4] + 2.*(ini.ftemp[1,:,3] + ini.ftemp[1,:,7] + ini.ftemp[1,:,6]))))
    ini.f[1,:,1] = (ini.ftemp[1,:,3] + (2./3.)*ini.ru[:])
    ini.f[1,:,5] = (ini.ftemp[1,:,7] + (1./6.)*ini.ru[:] - (1./2.)*(ini.ftemp[1,:,2] - ini.ftemp[1,:,4]))
    ini.f[1,:,8] = (ini.ftemp[1,:,6] + (1./6.)*ini.ru[:] - (1./2.)*(ini.ftemp[1,:,4] - ini.ftemp[1,:,2]))

    ini.f[1,:,1][ini.f[1,:,1] < 0] = 0
    ini.f[1,:,5][ini.f[1,:,5] < 0] = 0
    ini.f[1,:,8][ini.f[1,:,8] < 0] = 0

    
    # #Bottom left Corner
    # ini.f[1,1,1] = ini.f[1,-2,3]
    # ini.f[1,1,2] = ini.f[1,-2,4]
    # ini.f[1,1,5] = ini.f[1,-2,7]
    # ini.f[1,1,6] = 0.5*(ini.rho_left - (ini.f[1,-2,0] + 2*(ini.f[1,-2,1] + ini.f[1,-2,2] + ini.f[1,-2,5])))
    # ini.f[1,1,8] = ini.f[1,-2,6]

    # #Top left Corner
    # ini.f[1,-2,1] = ini.f[1,-2,3]
    # ini.f[1,-2,4] = ini.f[1,-2,2]
    # ini.f[1,-2,8] = ini.f[1,-2,6]
    # ini.f[1,-2,7] = 0.5*(ini.rho_left - (ini.f[1,-2,0] + 2*(ini.f[1,-2,1] + ini.f[1,-2,2] + ini.f[1,-2,6])))
    # ini.f[1,-2,5] = ini.f[1,-2,7]

    #Pressure BC on the right
    #OK
    ini.ru[:] = ini.rho0*(((ini.ftemp[-2,:,0] + ini.ftemp[-2,:,2] + ini.ftemp[-2,:,4] + 2.*(ini.ftemp[-2,:,1] + ini.ftemp[-2,:,5] + ini.ftemp[-2,:,8]))) - ini.rho_right)
    ini.f[-2,:,3] = (ini.ftemp[-2,:,1] - (2./3.)*ini.ru[:]) 
    ini.f[-2,:,7] = (ini.ftemp[-2,:,5] - (1./6.)*ini.ru[:] + (1./2.)*(ini.ftemp[-2,:,2] - ini.ftemp[-2,:,4]))
    ini.f[-2,:,6] = (ini.ftemp[-2,:,8] - (1./6.)*ini.ru[:] + (1./2.)*(ini.ftemp[-2,:,4] - ini.ftemp[-2,:,2]))

    ini.f[-2,:,3][ini.f[-2,:,3] < 0] = 0
    ini.f[-2,:,7][ini.f[-2,:,7] < 0] = 0
    ini.f[-2,:,6][ini.f[-2,:,6] < 0] = 0

    # print "ini.f[1,:,1] = ", ini.f[1,:,1]
    # print "ini.f[1,:,5] = ", ini.f[1,:,5]
    # print "ini.f[1,:,8] = ", ini.f[1,:,8]
    # print "ini.f[-2,:,3] = ", ini.f[-2,:,3]
    # print "ini.f[-2,:,7] = ", ini.f[-2,:,7]
    # print "ini.f[-2,:,6] = ", ini.f[-2,:,6]

    # ... computing density for imaging
    # ini.rho[:,:] = 0.
    # for a in range(9):
    #     ini.f[:,:,a] = np.where(abs(ini.f[:,:,a]) > ini.f_tol, ini.f[:,:,a], ini.f_tol)
    #     ini.rho[:,:] += np.where(ini.m[:,:] <> 1, ini.f[:,:,a], 0)
    # ini.rho[:,:] = ini.rho[:,:]*ini.rho0

    # print "ini.f[1,5,:] = ", ini.f[1,5,:]
    # print "ini.rho[2,5] = ", ini.rho[2,5]
    

    #Streaming step
    #ini.ftemp[:,:,:] = 0.

    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            for a in range(9):
                ini.ftemp[i,j,a] = 0
    
    alpha = 1.
    for j in range(1,ini.sizeY_+1):
        j_n = (j-1) #if j > 0 else (ini.sizeY_+1)
        j_p = (j+1) #if j < (ini.sizeY_ + 1) else 0 
        for i in range(1,ini.sizeX_+ 1):
            i_n = (i-1) 
            i_p = (i+1)
            if (ini.m[i,j] <> 1):
                #For streaming part, if adjacent grid is a m, then the density distribution will propagate, 
                #Else, the density is bounced back to the same grid, but with different direction        
                ini.ftemp[i,j,0] = ini.f[i,j,0]
                if (ini.m[i_p,j] <> 1):    ini.ftemp[i_p,j,1] = ini.f[i,j,1]
                else:                      ini.ftemp[i,j,3] = ini.f[i,j,1]

                if (ini.m[i,j_p] <> 1):    ini.ftemp[i,j_p,2] = ini.f[i,j,2]
                else:                      ini.ftemp[i,j,4] = ini.f[i,j,2]

                if (ini.m[i_n,j] <> 1):    ini.ftemp[i_n,j,3] = ini.f[i,j,3]
                else:                      ini.ftemp[i,j,1] = ini.f[i,j,3]

                if (ini.m[i,j_n] <> 1):    ini.ftemp[i,j_n,4] = ini.f[i,j,4]
                else:                      ini.ftemp[i,j,2] = ini.f[i,j,4]

                if (ini.m[i_p,j_p] <> 1):  ini.ftemp[i_p,j_p,5] = ini.f[i,j,5]
                else:
                    if (ini.m[i,j_p] == 1 and ini.m[i_p,j] == 1): ini.ftemp[i,j,7] += ini.f[i,j,5]
                    elif (ini.m[i,j_p] == 1 and ini.m[i_p,j] <> 1):
                        #Slip flow case 1 
                        if ((ini.rho[i_p,j]) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ini.rho[i_p,j]/ini.rho0)**ini.gamma)))
                        ini.ftemp[i,j,7] += alpha*ini.f[i,j,5]; ini.ftemp[i_p,j,8] += (1-alpha)*ini.f[i,j,5]
                    elif (ini.m[i,j_p] <> 1 and ini.m[i_p,j] == 1):
                        #Slip flow case 2
                        if ((ini.rho[i,j_p]) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ini.rho[i,j_p]/ini.rho0)**ini.gamma)))
                        ini.ftemp[i,j,7] += alpha*ini.f[i,j,5]; ini.ftemp[i,j_p,6] += (1-alpha)*ini.f[i,j,5]
                    elif (ini.m[i,j_p] <> 1 and ini.m[i_p,j] <> 1):
                        #Slip flow case 3
                        ave_rho = (ini.rho[i,j_p] + ini.rho[i_p,j])/2
                        if ((ave_rho) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ave_rho/ini.rho0)**ini.gamma)))
                        ini.ftemp[i,j,7] += alpha*ini.f[i,j,5]; ini.ftemp[i_p,j,8] += (1-alpha/2)*ini.f[i,j,5]; ini.ftemp[i,j_p,6] += (1-alpha/2)*ini.f[i,j,5]

                if (ini.m[i_n,j_p] <> 1):  ini.ftemp[i_n,j_p,6] = ini.f[i,j,6]
                else:
                    if (ini.m[i_n,j] == 1 and ini.m[i,j_p] == 1): ini.ftemp[i,j,8] += ini.f[i,j,6]
                    elif (ini.m[i_n,j] <> 1 and ini.m[i,j_p] == 1):
                        #Slip flow case 1
                        if ((ini.rho[i_n,j]) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ini.rho[i_n,j]/ini.rho0)**ini.gamma)))                      
                        ini.ftemp[i,j,8] += alpha*ini.f[i,j,6]; ini.ftemp[i_n,j,7] += (1-alpha)*ini.f[i,j,6]

                    elif (ini.m[i_n,j] == 1 and ini.m[i,j_p] <> 1):
                        #Slip flow case 2 
                        if ((ini.rho[i,j_p]) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ini.rho[i,j_p]/ini.rho0)**ini.gamma)))                      
                        ini.ftemp[i,j,8] += alpha*ini.f[i,j,6]; ini.ftemp[i,j_p,5] += (1-alpha)*ini.f[i,j,6]

                    elif (ini.m[i_n,j] <> 1 and ini.m[i,j_p] <> 1):
                        #Slip flow case 3
                        ave_rho = (ini.rho[i_n,j] + ini.rho[i,j_p])/2
                        if ((ave_rho) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ave_rho/ini.rho0)**ini.gamma)))
                        ini.ftemp[i,j,8] += alpha*ini.f[i,j,6]; ini.ftemp[i_n,j,7] += (1-alpha/2)*ini.f[i,j,6]; ini.ftemp[i,j_p,5] += (1-alpha/2)*ini.f[i,j,6]

                if (ini.m[i_n,j_n] <> 1):  ini.ftemp[i_n,j_n,7] = ini.f[i,j,7]
                else:
                    if (ini.m[i_n,j] == 1 and ini.m[i,j_n] == 1): ini.ftemp[i,j,5] += ini.f[i,j,7]
                    elif (ini.m[i_n,j] <> 1 and ini.m[i,j_n] == 1):
                        #Slip flow case 1
                        if ((ini.rho[i_n,j]) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ini.rho[i_n,j]/ini.rho0)**ini.gamma)))                      
                        ini.ftemp[i,j,5] += alpha*ini.f[i,j,7]; ini.ftemp[i_n,j,6] += (1-alpha)*ini.f[i,j,7]

                    elif (ini.m[i_n,j] == 1 and ini.m[i,j_n] <> 1):
                        #Slip flow case 2 
                        if ((ini.rho[i,j_n]) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ini.rho[i,j_n]/ini.rho0)**ini.gamma)))                      
                        ini.ftemp[i,j,5] += alpha*ini.f[i,j,7]; ini.ftemp[i,j_n,8] += (1-alpha)*ini.f[i,j,7]

                    elif (ini.m[i_n,j] <> 1 and ini.m[i,j_n] <> 1):
                        #Slip flow case 3
                        ave_rho = (ini.rho[i_n,j] + ini.rho[i,j_n])/2
                        if ((ave_rho) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ave_rho/ini.rho0)**ini.gamma)))
                        ini.ftemp[i,j,5] += alpha*ini.f[i,j,7]; ini.ftemp[i_n,j,6] += (1-alpha/2)*ini.f[i,j,7]; ini.ftemp[i,j_n,8] += (1-alpha/2)*ini.f[i,j,7]

                if (ini.m[i_p,j_n] <> 1):  ini.ftemp[i_p,j_n,8] = ini.f[i,j,8]
                else:
                    if (ini.m[i_p,j] == 1 and ini.m[i,j_n] == 1): ini.ftemp[i,j,6] += ini.f[i,j,8]
                    elif (ini.m[i_p,j] <> 1 and ini.m[i,j_n] == 1):
                        #Slip flow case 1
                        if ((ini.rho[i_p,j]) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ini.rho[i_p,j]/ini.rho0)**ini.gamma)))                      
                        ini.ftemp[i,j,6] += alpha*ini.f[i,j,8]; ini.ftemp[i_p,j,5] += (1-alpha)*ini.f[i,j,8]

                    elif (ini.m[i_p,j] == 1 and ini.m[i,j_n] <> 1):
                        #Slip flow case 2 
                        if ((ini.rho[i,j_n]) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ini.rho[i,j_n]/ini.rho0)**ini.gamma)))                      
                        ini.ftemp[i,j,6] += alpha*ini.f[i,j,8]; ini.ftemp[i,j_n,7] += (1-alpha)*ini.f[i,j,8]

                    elif (ini.m[i_p,j] <> 1 and ini.m[i,j_n] <> 1):
                        #Slip flow case 3
                        ave_rho = (ini.rho[i_p,j] + ini.rho[i,j_n])/2
                        if ((ave_rho) < 0.1): alpha = 0
                        else: alpha = 1 - (1/(ini.beta*((ave_rho/ini.rho0)**ini.gamma)))
                        ini.ftemp[i,j,6] += alpha*ini.f[i,j,8]; ini.ftemp[i_p,j,5] += (1-alpha/2)*ini.f[i,j,8]; ini.ftemp[i,j_n,7] += (1-alpha/2)*ini.f[i,j,8]

    print "alpha = ", alpha
    # ... and then computing macroscopic density and velocity for each lattice point, after shifting
    ini.rho[:,:] = 0.
    ini.ux[:,:] = 0.
    ini.uy[:,:] = 0.
    for a in range(9):
        ini.rho[:,:] += ini.ftemp[:,:,a]    
        ini.ux[:,:] += ini.e_[0,a]*ini.ftemp[:,:,a]
        ini.uy[:,:] += ini.e_[1,a]*ini.ftemp[:,:,a]
    ini.rho[:,:] = ini.rho[:,:]*ini.rho0
    ini.ux[:,:] = ini.ux[:,:]/ini.rho[:,:]
    ini.uy[:,:] = ini.uy[:,:]/ini.rho[:,:]
    ini.u[:,:] = sqrt((ini.ux[:,:]**2 + ini.uy[:,:]**2)/2)

    # if mod(t,100) == 0:
    #     plt.figure(2)
    #     # plt.plot(ini.rho[1:-1,5])
    #     # plt.plot(ini.ux[-3,:],label='Outlet')
    #     plt.plot(ini.ux[(ini.sizeX_+2)//2,:],label='Middle')
    #     plt.ylim(bottom = 0, top = 0.5)
    #     plt.xlim(left = 0, right = 14)
    #     # plt.plot(ini.ux[2,:],label='Inlet')
    #     plt.legend(loc='upper right')
    #     plt.show
    #     plt.pause(0.001)
    #     plt.clf()

    # print ini.ux[(ini.sizeX_+2)//2,:]


    #Computing equilibrium distribution function
    ini.fct1[:,:] = ini.w[0]*ini.rho[:,:]/ini.rho0
    ini.fct2[:,:] = ini.w[1]*ini.rho[:,:]/ini.rho0
    ini.fct3[:,:] = ini.w[2]*ini.rho[:,:]/ini.rho0

    ini.uxeq[:,:] = ini.ux[:,:]                                       #uxeq will incorporate external forces, if any
    ini.uyeq[:,:] = ini.uy[:,:] 

    ini.uxsq[:,:] = ini.uxeq[:,:]*ini.uxeq[:,:]
    ini.uysq[:,:] = ini.uyeq[:,:]*ini.uyeq[:,:]

    ini.uxuy5[:,:] = ini.uxeq[:,:] + ini.uyeq[:,:]
    ini.uxuy6[:,:] = -ini.uxeq[:,:] + ini.uyeq[:,:]
    ini.uxuy7[:,:] = -ini.uxeq[:,:] - ini.uyeq[:,:]
    ini.uxuy8[:,:] = ini.uxeq[:,:] - ini.uyeq[:,:]

    ini.usq[:,:] = ini.uxsq[:,:] + ini.uysq[:,:]

    ini.feq[:,:,0] = ini.fct1[:,:]*(1.                                                   - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,1] = ini.fct2[:,:]*(1. + ini.c_eq[0]*ini.uxeq[:,:] + ini.c_eq[1]*ini.uxsq[:,:]    - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,2] = ini.fct2[:,:]*(1. + ini.c_eq[0]*ini.uyeq[:,:] + ini.c_eq[1]*ini.uysq[:,:]    - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,3] = ini.fct2[:,:]*(1. - ini.c_eq[0]*ini.uxeq[:,:] + ini.c_eq[1]*ini.uxsq[:,:]    - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,4] = ini.fct2[:,:]*(1. - ini.c_eq[0]*ini.uyeq[:,:] + ini.c_eq[1]*ini.uysq[:,:]    - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,5] = ini.fct3[:,:]*(1. + ini.c_eq[0]*ini.uxuy5[:,:] + ini.c_eq[1]*ini.uxuy5[:,:]*ini.uxuy5[:,:]  - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,6] = ini.fct3[:,:]*(1. + ini.c_eq[0]*ini.uxuy6[:,:] + ini.c_eq[1]*ini.uxuy6[:,:]*ini.uxuy6[:,:]  - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,7] = ini.fct3[:,:]*(1. + ini.c_eq[0]*ini.uxuy7[:,:] + ini.c_eq[1]*ini.uxuy7[:,:]*ini.uxuy7[:,:]  - ini.c_eq[2]*ini.usq[:,:])
    ini.feq[:,:,8] = ini.fct3[:,:]*(1. + ini.c_eq[0]*ini.uxuy8[:,:] + ini.c_eq[1]*ini.uxuy8[:,:]*ini.uxuy8[:,:]  - ini.c_eq[2]*ini.usq[:,:])

    # print "ini.ftemp[2,5,:] = ",ini.ftemp[2,5,:]
    # print "ini.feq[2,5,:] = ",ini.feq[2,5,:]

    #Collision step
    #ONLY TEST BEFORE FRIDAY PRAYER 19/10
    ini.tau[:,:,:] = maximum(ini.tau0, 1 - (ini.feq[:,:,:]/ini.ftemp[:,:,:]))
    for a in range(9):
       ini.f[:,:,a] = np.where(ini.m[:,:] == 0, ini.ftemp[:,:,a] - (ini.ftemp[:,:,a] - ini.feq[:,:,a]) / ini.tau[:,:,a], ini.f[:,:,a])

    #print "ini.f[2,5,:] = ", ini.f[2,5,:]

    #Plotting 
    print "Time = ", t
    print "Mass = ", sum(ini.rho)
    print "Velocity x dir = ", sum(ini.ux)
    print "Velocity y dir = ", sum(ini.uy)
    # print "rho_inlet = ", ini.rho[2,(ini.sizeY_+2)//2]
    # print "vel_inlet = ", ini.ux[2,(ini.sizeY_+2)//2]
    # print "rho_outlet = ", ini.rho[-3,(ini.sizeY_+2)//2]
    # print "vel_outlet = ", ini.ux[-3,(ini.sizeY_+2)//2]
    # print "rho_up = ", ini.rho[(ini.sizeX_+2)//2 - 2,(ini.sizeY_+2)//2]
    # print "rho_down = ", ini.rho[(ini.sizeX_+2)//2 + 2,(ini.sizeY_+2)//2]
    # print "vel_up = ", ini.ux[(ini.sizeX_+2)//2 - 2,(ini.sizeY_+2)//2]
    # print "vel_down = ", ini.ux[(ini.sizeX_+2)//2 + 2,(ini.sizeY_+2)//2]
    #print ini.f[2:ini.sizeX_,2:ini.sizeY_,:].min()
    #print ini.ux[ini.sizeX_//2,ini.sizeY_//2] 
    print "Middle throughput = ", sum(ini.ux[(ini.sizeX_+2)//2,2:-3]*ini.rho[(ini.sizeX_+2)//2,2:-3])/(ini.sizeY_-2)

    #To save the density and velocity distribution
    # np.save(os.path.join(ini.name, "rho_" + ini.name + "_" + str(t).zfill(4)), ini.rho)
    # np.save(os.path.join(ini.name, "ux_" + ini.name + "_" + str(t).zfill(4)), ini.ux)     
    
    if mod(t,50) == 0:
        varm = ini.rho.transpose()        #Change the variable to the one that will be plotted: rho, ux, or uy
        plt.figure(1)
        ax = sns.heatmap(varm, annot=False, vmin=0., vmax=1.5, cmap='RdYlBu_r', mask=ini.m.transpose())           #For ux
        #ax = sns.heatmap(varm, annot=False, vmin=-0.5, vmax=0.5, cmap='RdYlBu_r', mask=ini.m.transpose())           #For rho
        #figure(num=1, figsize=(15,6), dpi=80, facecolor='w', edgecolor='k')
        fig = plt.gcf()
        fig.set_size_inches(10.,5., forward=True)
        #ax = sns.heatmap(varm, annot=False, vmin=5.*ini.f_init, vmax=18*ini.f_init, cmap='RdYlBu_r')           #For rho
        ax.invert_yaxis()
        plt.pause(0.001)
        #plt.savefig("vel_"+ini.name2+"_"+str(t).zfill(6)+".png")
        plt.clf()

    #Plot cross-section velocity


    
    #Plot drag force vs time
    # if mod(t,5) == 0:
    #     plt.figure(3)
    #     plt.plot(range(ini.T), ini.Ft, label='Outlet')
    #     plt.legend(loc='upper right')
    #     plt.show
    #     plt.pause(0.001)
    #     plt.clf()

#To save the force
# np.save(os.path.join(ini.name, "dragForce_" + ini.name), ini.Ft)

####################################################### OUTPUT ###########################################################################