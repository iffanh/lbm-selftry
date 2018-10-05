#This is an attempt to make a program to simulate fluid flow using LBM
#This is for single-phase flow

####################################################### LIBRARIES ###########################################################################
import math
import numpy as np
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import lbm_input as ini             #import from theinput file

####################################################### SIMULATION ###########################################################################

for t in range(ini.T):

    #Zou and He velocity BCs on west side
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.m[i][j] == 2:
                ini.rho[i][j] = (ini.f[i][j][0] + ini.f[i][j][2] + ini.f[i][j][4] + 2.*(ini.f[i][j][3] + ini.f[i][j][7] + ini.f[i][j][6])) / (1 - ini.ux0)
                ru = ini.rho[i][j]*ini.ux0
                ini.f[i][j][1] = ini.f[i][j][3] + (2./3.)*ru
                ini.f[i][j][5] = ini.f[i][j][7] + (1./6.)*ru - (1./2.)*(ini.f[i][j][2] - ini.f[i][j][4])
                ini.f[i][j][8] = ini.f[i][j][6] + (1./6.)*ru - (1./2.)*(ini.f[i][j][4] - ini.f[i][j][2])
                

    # ... computing density for imaging        
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.m[i][j] == 0:
                ini.rho[i][j] = 0
                for a in range(9):
                    ini.ftemp[i][j][a] = 0
                    if ini.f[i][j][a] < ini.f_tol: ini.f[i][j][a] = 0
                    ini.rho[i][j] += ini.f[i][j][a]    

    #Streaming step
    for j in range(1,ini.sizeY_+ 1):
        j_n = (j-1) 
        j_p = (j+1)
        for i in range(1,ini.sizeX_+ 1):
            i_n = (i-1) 
            i_p = (i+1)
            if (ini.m[i][j] <> 1):
                #For streaming part, if adjacent grid is a m, then the density distribution will propagate, 
                #Else, the density is bounced back to the same grid, but with different direction
                
                ini.ftemp[i][j][0] = ini.f[i][j][0]
                
                if (ini.m[i_p][j] <> 1):    ini.ftemp[i_p][j][1] = ini.f[i][j][1]#; f[i][j][1] = 0
                else:                           ini.ftemp[i][j][3] = ini.f[i][j][1]#; f[i][j][1] = 0

                if (ini.m[i][j_p] <> 1):    ini.ftemp[i][j_p][2] = ini.f[i][j][2]#; f[i][j][2] = 0
                else:                           ini.ftemp[i][j][4] = ini.f[i][j][2]#; f[i][j][2] = 0

                if (ini.m[i_n][j] <> 1):    ini.ftemp[i_n][j][3] = ini.f[i][j][3]#; f[i][j][3] = 0 
                else:                           ini.ftemp[i][j][1] = ini.f[i][j][3]#; f[i][j][3] = 0 

                if (ini.m[i][j_n] <> 1):    ini.ftemp[i][j_n][4] = ini.f[i][j][4]#; f[i][j][4] = 0
                else:                           ini.ftemp[i][j][2] = ini.f[i][j][4]#; f[i][j][4] = 0

                if (ini.m[i_p][j_p] <> 1):  ini.ftemp[i_p][j_p][5] = ini.f[i][j][5]#; f[i][j][5] = 0
                else:                           ini.ftemp[i][j][7] = ini.f[i][j][5]#; f[i][j][5] = 0

                if (ini.m[i_n][j_p] <> 1):  ini.ftemp[i_n][j_p][6] = ini.f[i][j][6]#; f[i][j][6] = 0
                else:                           ini.ftemp[i][j][8] = ini.f[i][j][6]#; f[i][j][6] = 0

                if (ini.m[i_n][j_n] <> 1):  ini.ftemp[i_n][j_n][7] = ini.f[i][j][7]#; f[i][j][7] = 0
                else:                           ini.ftemp[i][j][5] = ini.f[i][j][7]#; f[i][j][7] = 0

                if (ini.m[i_p][j_n] <> 1):  ini.ftemp[i_p][j_n][8] = ini.f[i][j][8]#; f[i][j][8] = 0
                else:                           ini.ftemp[i][j][6] = ini.f[i][j][8]#; f[i][j][8] = 0
    
    # ... and then computing macroscopic density and velocity for each lattice point, after shifting        
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.m[i][j] <> 1:# or m[i][j] == 1:
                ini.rho[i][j] = 0
                ini.ux[i][j] = 0
                ini.uy[i][j] = 0
                for a in range(9):
                    ini.rho[i][j] += ini.f[i][j][a]    
                    ini.ux[i][j] += ini.e_x[a]*ini.f[i][j][a]
                    ini.uy[i][j] += ini.e_y[a]*ini.f[i][j][a]
                ini.ux[i][j] = ini.ux[i][j]/ini.rho[i][j] if ini.rho[i][j] <> 0 else 0
                ini.uy[i][j] = ini.uy[i][j]/ini.rho[i][j] if ini.rho[i][j] <> 0 else 0

    #Computing equilibrium distribution function
    for j in range(1,ini.sizeY_+ 1):
        for i in range(1,ini.sizeX_+ 1):
            if ini.m[i][j] == 0:
                fct1 = ini.w1*ini.rho[i][j]
                fct2 = ini.w2*ini.rho[i][j]
                fct3 = ini.w3*ini.rho[i][j]

                ini.uxeq[i][j] = ini.ux[i][j]                                       #uxeq will incorporate external forces, if any
                ini.uyeq[i][j] = ini.uy[i][j] 

                uxsq = ini.uxeq[i][j]*ini.uxeq[i][j]
                uysq = ini.uyeq[i][j]*ini.uyeq[i][j]

                uxuy5 = ini.uxeq[i][j] + ini.uyeq[i][j]
                uxuy6 = -ini.uxeq[i][j] + ini.uyeq[i][j]
                uxuy7 = -ini.uxeq[i][j] - ini.uyeq[i][j]
                uxuy8 = ini.uxeq[i][j] - ini.uyeq[i][j]

                usq = uxsq + uysq

                ini.feq[i][j][0] = fct1*(1.                              - ini.f3*usq)
                ini.feq[i][j][1] = fct2*(1. + ini.f1*ini.uxeq[i][j] + ini.f2*uxsq    - ini.f3*usq)
                ini.feq[i][j][2] = fct2*(1. + ini.f1*ini.uyeq[i][j] + ini.f2*uysq    - ini.f3*usq)
                ini.feq[i][j][3] = fct2*(1. - ini.f1*ini.uxeq[i][j] + ini.f2*uxsq    - ini.f3*usq)
                ini.feq[i][j][4] = fct2*(1. - ini.f1*ini.uyeq[i][j] + ini.f2*uysq    - ini.f3*usq)
                ini.feq[i][j][5] = fct3*(1. + ini.f1*uxuy5 + ini.f2*uxuy5*uxuy5  - ini.f3*usq)
                ini.feq[i][j][6] = fct3*(1. + ini.f1*uxuy6 + ini.f2*uxuy6*uxuy6  - ini.f3*usq)
                ini.feq[i][j][7] = fct3*(1. + ini.f1*uxuy7 + ini.f2*uxuy7*uxuy7  - ini.f3*usq)
                ini.feq[i][j][8] = fct3*(1. + ini.f1*uxuy8 + ini.f2*uxuy8*uxuy8  - ini.f3*usq)

    #Collision step
    for j in range(1,ini.sizeY_+1):
        j_n = (j-1) 
        j_p = (j+1)
        for i in range(1,ini.sizeX_+1):
            i_n = (i-1) 
            i_p = (i+1)
            if ini.m[i][j] == 0:
                for a in range(9):
                    #Only calculate the one with no bounce back? Not yet implemented
                    ini.f[i][j][a] = ini.ftemp[i][j][a] - (ini.ftemp[i][j][a] - ini.feq[i][j][a]) / ini.tau

    densityM = zip(*ini.rho)         #Transpose matrix rho
    print "Mass = ", np.sum(densityM)
    print "Velocity x dir = ", np.sum(ini.ux)
    print "Velocity y dir = ", np.sum(ini.uy)
    ax = sns.heatmap(densityM, annot=False, vmin=0, vmax=5)
    ax.invert_yaxis()
    plt.draw()
    plt.pause(0.1)
    plt.clf()
####################################################### OUTPUT ###########################################################################

        
