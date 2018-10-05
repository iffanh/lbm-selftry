#This is an attempt to make a program to simulate fluid flow using LBM
#This is for single-phase flow

####################################################### LIBRARIES ###########################################################################
import math
import numpy as np
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import lbm_input as ini

# ####################################################### PREPARATION ###########################################################################
# #The size of grid
# sizeX_ = 40         #length in x-direction
# sizeY_ = 16         #length in y-direction

# #The number of iteration
# T = 200            #Total time used in the simulation
# dt = 1             #time interval

# solid = [[0 for j in xrange(sizeY_ + 2)] for i in xrange(sizeX_ + 2)]                       #Presence of solid or not, 1 means solid

# #Declaring variables

# rho = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]                         #Density of the lattice point, 
# rho_temp = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]
# ux = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]                          #Macroscopic velocity of the lattice point 
# uy = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]
# uxeq = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]                        #Macroscopic velocity of the lattice point 
# uyeq = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]
# f = [[[0 for k in xrange(9)] for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]      #Density distribution of the a point f[x position][y position][index]
# ftemp = [[[0 for k in xrange(9)] for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]
# feq = [[[0 for k in xrange(9)] for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]

# #Constants used
# tau = 1.5
# e_x = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0]          
# e_y = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]
# w = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]
# w1 = 4./9.
# w2 = 1./9.
# w3 = 1./36.
# f1 = 3.
# f2 = 9./2.
# f3 = 3./2.
# f_tol = 0. 

# ####################################################### FUNCTIONS ###########################################################################
# #To find out whether a solid is in contact in fluid or not
# def is_interior(i,j):       
#     if solid[i][j] == 1 and solid[i-1][j] == 1 and solid[i+1][j] == 1 and solid[i][j-1] == 1 and solid[i][j+1] == 1 and solid[i-1][j-1] == 1 and solid[i+1][j-1] == 1 and solid[i-1][j+1] == 1 and solid[i+1][j+1] == 1:
#         return True

# ####################################################### INTIALIZATION ###########################################################################

# ####Sacrificed grid 
# for i in range(0,sizeX_+ 2):
#     for j in range(0,sizeY_ + 2):
#         solid[i][0] = 0
#         solid[i][sizeY_+1] = 0
#         solid[0][j] = 0
#         solid[sizeX_+1][j] = 0

# ####Solid presence in Grid 
# #The outer boundary
# #0 is a normal no solid boundary
# #1 is a bounce back boundary
# #2 is a von neumann boundary
# for i in range(1,sizeX_+ 1):
#     for j in range(1,sizeY_ + 1):
#         solid[i][1] = 1
#         solid[i][sizeY_] = 1
#         solid[1][j] = 2
#         solid[sizeX_][j] = 0

# solid[1][1] = 1
# solid[1][sizeY_] = 1
# solid[sizeX_][1] = 1
# solid[sizeX_][sizeY_] = 1

# #solid[sizeX_][2] = 1
# #solid[sizeX_][sizeY_-1] = 1

# #The cylinder
# solid[7][10] = 1; solid[8][10] = 1; solid[9][10] = 1
# solid[6][9] = 1; solid[7][9] = 1; solid[8][9] = 1; solid[9][9] = 1; solid[10][9] = 1
# solid[6][8] = 1; solid[7][8] = 1; solid[8][8] = 1; solid[9][8] = 1; solid[10][8] = 1
# solid[6][7] = 1; solid[7][7] = 1; solid[8][7] = 1; solid[9][7] = 1; solid[10][7] = 1
# solid[7][6] = 1; solid[8][6] = 1; solid[9][6] = 1

# ####Initial Condition for density distribution, f
# #Initialize density distribution f, ...

# f_init = 0.1
# for j in range(1,sizeY_+ 1):
#     for i in range(1, sizeX_+ 1):                 
#         if solid[i][j] == 0:
#             for a in range(9):
#                 f[i][j][a] = f_init 

# ###Von Neumann Boundary condition
# #Initializing flux boundary density distribution
# ux0 = 0.8
# for i in range(1,sizeX_+ 1):
#     for j in range(1,sizeY_ + 1):
#         if solid[i][j] == 2:
#             #West side
#             for a in [0,2,3,4,6,7]:
#                 f[i][j][a] = f_init

# #f[10][10][5] = 10.


####################################################### SIMULATION ###########################################################################

for t in range(ini.T):

    #Zou and He velocity BCs on west side
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.solid[i][j] == 2:
                ini.rho[i][j] = (ini.f[i][j][0] + ini.f[i][j][2] + ini.f[i][j][4] + 2.*(ini.f[i][j][3] + ini.f[i][j][7] + ini.f[i][j][6])) / (1 - ini.ux0)
                ru = ini.rho[i][j]*ini.ux0
                ini.f[i][j][1] = ini.f[i][j][3] + (2./3.)*ru
                ini.f[i][j][5] = ini.f[i][j][7] + (1./6.)*ru - (1./2.)*(ini.f[i][j][2] - ini.f[i][j][4])
                ini.f[i][j][8] = ini.f[i][j][6] + (1./6.)*ru - (1./2.)*(ini.f[i][j][4] - ini.f[i][j][2])
                

    # ... computing density for imaging        
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.solid[i][j] == 0:
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
            if (ini.solid[i][j] <> 1):
                #For streaming part, if adjacent grid is a solid, then the density distribution will propagate, 
                #Else, the density is bounced back to the same grid, but with different direction
                
                ini.ftemp[i][j][0] = ini.f[i][j][0]
                
                if (ini.solid[i_p][j] <> 1):    ini.ftemp[i_p][j][1] = ini.f[i][j][1]#; f[i][j][1] = 0
                else:                           ini.ftemp[i][j][3] = ini.f[i][j][1]#; f[i][j][1] = 0

                if (ini.solid[i][j_p] <> 1):    ini.ftemp[i][j_p][2] = ini.f[i][j][2]#; f[i][j][2] = 0
                else:                           ini.ftemp[i][j][4] = ini.f[i][j][2]#; f[i][j][2] = 0

                if (ini.solid[i_n][j] <> 1):    ini.ftemp[i_n][j][3] = ini.f[i][j][3]#; f[i][j][3] = 0 
                else:                           ini.ftemp[i][j][1] = ini.f[i][j][3]#; f[i][j][3] = 0 

                if (ini.solid[i][j_n] <> 1):    ini.ftemp[i][j_n][4] = ini.f[i][j][4]#; f[i][j][4] = 0
                else:                           ini.ftemp[i][j][2] = ini.f[i][j][4]#; f[i][j][4] = 0

                if (ini.solid[i_p][j_p] <> 1):  ini.ftemp[i_p][j_p][5] = ini.f[i][j][5]#; f[i][j][5] = 0
                else:                           ini.ftemp[i][j][7] = ini.f[i][j][5]#; f[i][j][5] = 0

                if (ini.solid[i_n][j_p] <> 1):  ini.ftemp[i_n][j_p][6] = ini.f[i][j][6]#; f[i][j][6] = 0
                else:                           ini.ftemp[i][j][8] = ini.f[i][j][6]#; f[i][j][6] = 0

                if (ini.solid[i_n][j_n] <> 1):  ini.ftemp[i_n][j_n][7] = ini.f[i][j][7]#; f[i][j][7] = 0
                else:                           ini.ftemp[i][j][5] = ini.f[i][j][7]#; f[i][j][7] = 0

                if (ini.solid[i_p][j_n] <> 1):  ini.ftemp[i_p][j_n][8] = ini.f[i][j][8]#; f[i][j][8] = 0
                else:                           ini.ftemp[i][j][6] = ini.f[i][j][8]#; f[i][j][8] = 0
    
    # ... and then computing macroscopic density and velocity for each lattice point, after shifting        
    for j in range(1,ini.sizeY_+1):
        for i in range(1,ini.sizeX_+1):
            if ini.solid[i][j] <> 1:# or solid[i][j] == 1:
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
            if ini.solid[i][j] == 0:
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
            if ini.solid[i][j] == 0:
                for a in range(9):
                    #Only calculate the one with no bounce back? Not yet implemented
                    ini.f[i][j][a] = ini.ftemp[i][j][a] - (ini.ftemp[i][j][a] - ini.feq[i][j][a]) / ini.tau

                # f[i][j][0] = ftemp[i][j][0] - (ftemp[i][j][0] - feq[i][j][0]) / tau    
                # if (solid[i_p][j] <> 1):    f[i][j][1] = ftemp[i][j][1] - (ftemp[i][j][1] - feq[i][j][1]) / tau
                # else:                       f[i][j][3] = ftemp[i][j][3]
                # if (solid[i][j_p] <> 1):    f[i][j][2] = ftemp[i][j][2] - (ftemp[i][j][2] - feq[i][j][2]) / tau
                # else:                       f[i][j][4] = ftemp[i][j][4]
                # if (solid[i_n][j] <> 1):    f[i][j][3] = ftemp[i][j][3] - (ftemp[i][j][3] - feq[i][j][3]) / tau
                # else:                       f[i][j][1] = ftemp[i][j][1]
                # if (solid[i][j_n] <> 1):    f[i][j][4] = ftemp[i][j][4] - (ftemp[i][j][4] - feq[i][j][4]) / tau
                # else:                       f[i][j][2] = ftemp[i][j][2]
                # if (solid[i_p][j_p] <> 1):  f[i][j][5] = ftemp[i][j][5] - (ftemp[i][j][5] - feq[i][j][5]) / tau
                # else:                       f[i][j][7] = ftemp[i][j][7]
                # if (solid[i_n][j_p] <> 1):  f[i][j][6] = ftemp[i][j][6] - (ftemp[i][j][6] - feq[i][j][6]) / tau
                # else:                       f[i][j][8] = ftemp[i][j][8]
                # if (solid[i_n][j_n] <> 1):  f[i][j][7] = ftemp[i][j][7] - (ftemp[i][j][7] - feq[i][j][7]) / tau
                # else:                       f[i][j][5] = ftemp[i][j][5]
                # if (solid[i_p][j_n] <> 1):  f[i][j][8] = ftemp[i][j][8] - (ftemp[i][j][8] - feq[i][j][8]) / tau 
                # else:                       f[i][j][6] = ftemp[i][j][6]

    densityM = zip(*ini.rho)         #Transpose matrix rho
    #densityM = zip(*solid)         #Transpose matrix rho
    print "Mass = ", np.sum(densityM)
    print "Velocity x dir = ", np.sum(ini.ux)
    print "Velocity y dir = ", np.sum(ini.uy)
    ax = sns.heatmap(densityM, annot=False, vmin=0, vmax=5)
    ax.invert_yaxis()
    plt.draw()
    plt.pause(0.1)
    plt.clf()
####################################################### OUTPUT ###########################################################################

        
