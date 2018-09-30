#This is an attempt to make a program to simulate fluid flow using LBM
#This is for single-phase flow

####################################################### LIBRARIES ###########################################################################
import math
import numpy as np
import seaborn as sns; sns.set()
import matplotlib.pyplot as plt

####################################################### PREPARATION ###########################################################################
#The size of grid
sizeX_ = 40         #length in x-direction
sizeY_ = 16         #length in y-direction

#The number of iteration
T = 200            #Total time used in the simulation
dt = 1             #time interval

solid = [[0 for j in xrange(sizeY_ + 2)] for i in xrange(sizeX_ + 2)]                       #Presence of solid or not, 1 means solid

#Declaring variables

rho = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]                         #Density of the lattice point, 
rho_temp = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]
ux = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]                          #Macroscopic velocity of the lattice point 
uy = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]
uxeq = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]                        #Macroscopic velocity of the lattice point 
uyeq = [[0 for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]
f = [[[0 for k in xrange(9)] for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]      #Density distribution of the a point f[x position][y position][index]
ftemp = [[[0 for k in xrange(9)] for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]
feq = [[[0 for k in xrange(9)] for j in xrange(sizeY_+ 2)] for i in xrange(sizeX_+ 2)]

#Constants used
tau = 12.
e_x = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0]          
e_y = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]
w = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]
w1 = 4./9.
w2 = 1./9.
w3 = 1./36.
f1 = 3.
f2 = 9./2.
f3 = 3./2. 

####################################################### FUNCTIONS ###########################################################################
#To find out whether a solid is in contact in fluid or not
def is_interior(i,j):       
    if solid[i][j] == 1 and solid[i-1][j] == 1 and solid[i+1][j] == 1 and solid[i][j-1] == 1 and solid[i][j+1] == 1 and solid[i-1][j-1] == 1 and solid[i+1][j-1] == 1 and solid[i-1][j+1] == 1 and solid[i+1][j+1] == 1:
        return True

####################################################### INTIALIZATION ###########################################################################

####Solid presence in Grid 

for i in range(0,sizeX_+ 2):
    for j in range(0,sizeY_ + 2):
        solid[i][0] = 1
        solid[i][sizeY_+1] = 1
        solid[0][j] = 1
        solid[sizeX_+1][j] = 1

# for i in range(3*sizeX_//7,4*sizeX_//7):
#     for j in range(3*sizeY_//7,4*sizeY_//7):
#         solid[i][j] = 1

solid[7][10] = 1; solid[8][10] = 1; solid[9][10] = 1
solid[6][9] = 1; solid[7][9] = 1; solid[8][9] = 1; solid[9][9] = 1; solid[10][9] = 1
solid[6][8] = 1; solid[7][8] = 1; solid[8][8] = 1; solid[9][8] = 1; solid[10][8] = 1
solid[6][7] = 1; solid[7][7] = 1; solid[8][7] = 1; solid[9][7] = 1; solid[10][7] = 1
solid[7][6] = 1; solid[8][6] = 1; solid[9][6] = 1

####Initial Condition for density distribution, f
#Initialize density distribution f, ...

f_init = 0.2
for j in range(1,sizeY_+ 1):
    for i in range(1, sizeX_+ 1):                 
        if solid[i][j] == 0: #or solid[i][j] == 1:
            #f[i][j][0] = f_init  
            for a in range(9):
                f[i][j][a] = f_init #if i < ((sizeX_+2)//2) else f_init
        # #Initialize bottom left distribution
        # elif solid[i][j] == 1 and i == 1 and j == 1: f[i][j][5] = f_init                                
        # elif solid[i][j] == 1 and i == 2 and j == 1: f[i][j][5] = f_init; f[i][j][2] = f_init
        # elif solid[i][j] == 1 and i == 1 and j == 2: f[i][j][5] = f_init; f[i][j][1] = f_init
        # #Initialize upper left distribution
        # elif solid[i][j] == 1 and i == 1 and j == sizeY_: f[i][j][8] = f_init                           
        # elif solid[i][j] == 1 and i == 2 and j == sizeY_: f[i][j][8] = f_init; f[i][j][4] = f_init
        # elif solid[i][j] == 1 and i == 1 and j == sizeY_-1: f[i][j][8] = f_init; f[i][j][1] = f_init
        # #Initialize bottom right distribution
        # elif solid[i][j] == 1 and i == sizeX_ and j == 1: f[i][j][6] = f_init                           
        # elif solid[i][j] == 1 and i == sizeX_-1 and j == 1: f[i][j][6] = f_init; f[i][j][2] = f_init
        # elif solid[i][j] == 1 and i == sizeX_ and j == 2: f[i][j][6] = f_init; f[i][j][3] = f_init
        # #Initialize upper right distribution
        # elif solid[i][j] == 1 and i == sizeX_ and j == sizeY_: f[i][j][7] = f_init                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
        # elif solid[i][j] == 1 and i == sizeX_-1 and j == sizeY_: f[i][j][7] = f_init; f[i][j][4] = f_init
        # elif solid[i][j] == 1 and i == sizeX_ and j == sizeY_-1: f[i][j][7] = f_init; f[i][j][3] = f_init
        # #These four parts are for the side boundaries
        # elif solid[i][j] == 1 and i == 1: f[i][j][1] = f_init; f[i][j][5] = f_init; f[i][j][8] = f_init 
        # elif solid[i][j] == 1 and i == sizeX_: f[i][j][3] = f_init; f[i][j][6] = f_init; f[i][j][7] = f_init
        # elif solid[i][j] == 1 and j == 1: f[i][j][2] = f_init; f[i][j][5] = f_init; f[i][j][6] = f_init
        # elif solid[i][j] == 1 and j == sizeY_: f[i][j][4] = f_init; f[i][j][7] = f_init; f[i][j][8] = f_init
            


#f[1][5][1] = 1.0
f[20][7][3] = 10.

####################################################### SIMULATION ###########################################################################

for t in range(T):
    # ... and then computing macroscopic density and velocity for each lattice point        
    for j in range(1,sizeY_+1):
        for i in range(1,sizeX_+1):
            if solid[i][j] == 0:# or solid[i][j] == 1:
                rho[i][j] = 0
                #ux[i][j] = 0
                #uy[i][j] = 0
                for a in range(9):
                    rho[i][j] += f[i][j][a]    
                    #ux[i][j] += e_x[a]*f[i][j][a]
                    #uy[i][j] += e_y[a]*f[i][j][a]
                #ux[i][j] = ux[i][j]/rho[i][j] if rho[i][j] <> 0 else 0
                #uy[i][j] = ux[i][j]/rho[i][j] if rho[i][j] <> 0 else 0

    #Streaming step
    for j in range(1,sizeY_+ 1):
        j_n = (j-1) 
        j_p = (j+1)
        for i in range(1,sizeX_+ 1):
            i_n = (i-1) 
            i_p = (i+1)
            if (solid[i][j] == 0):# or (solid[i][j] == 1 and ((i == 1) or (i == sizeX_) or (j == 1) or (j == sizeY_)))):
                if (solid[i][j] == 0):      ftemp[i][j][0] = f[i][j][0]

                if (solid[i_p][j] == 0):    ftemp[i_p][j][1] = f[i][j][1]#; f[i][j][1] = 0
                else:                       ftemp[i][j][3] = f[i][j][1]#; f[i][j][1] = 0

                if (solid[i][j_p] == 0):    ftemp[i][j_p][2] = f[i][j][2]#; f[i][j][2] = 0
                else:                       ftemp[i][j][4] = f[i][j][2]#; f[i][j][2] = 0

                if (solid[i_n][j] == 0):    ftemp[i_n][j][3] = f[i][j][3]#; f[i][j][3] = 0 
                else:                       ftemp[i][j][1] = f[i][j][3]#; f[i][j][3] = 0 

                if (solid[i][j_n] == 0):    ftemp[i][j_n][4] = f[i][j][4]#; f[i][j][4] = 0
                else:                       ftemp[i][j][2] = f[i][j][4]#; f[i][j][4] = 0

                if (solid[i_p][j_p] == 0):  ftemp[i_p][j_p][5] = f[i][j][5]#; f[i][j][5] = 0
                else:                       ftemp[i][j][7] = f[i][j][5]#; f[i][j][5] = 0

                if (solid[i_n][j_p] == 0):  ftemp[i_n][j_p][6] = f[i][j][6]#; f[i][j][6] = 0
                else:                       ftemp[i][j][8] = f[i][j][6]#; f[i][j][6] = 0

                if (solid[i_n][j_n] == 0):  ftemp[i_n][j_n][7] = f[i][j][7]#; f[i][j][7] = 0
                else:                       ftemp[i][j][5] = f[i][j][7]#; f[i][j][7] = 0

                if (solid[i_p][j_n] == 0):  ftemp[i_p][j_n][8] = f[i][j][8]#; f[i][j][8] = 0
                else:                       ftemp[i][j][6] = f[i][j][8]#; f[i][j][8] = 0

                # ftemp[i][j][0] = f[i][j][0]
                # ftemp[i_p][j][1] = f[i][j][1]#; f[i][j][1] = 0
                # ftemp[i][j_p][2] = f[i][j][2]#; f[i][j][2] = 0
                # ftemp[i_n][j][3] = f[i][j][3]#; f[i][j][3] = 0 
                # ftemp[i][j_n][4] = f[i][j][4]#; f[i][j][4] = 0
                # ftemp[i_p][j_p][5] = f[i][j][5]#; f[i][j][5] = 0
                # ftemp[i_n][j_p][6] = f[i][j][6]#; f[i][j][6] = 0
                # ftemp[i_n][j_n][7] = f[i][j][7]#; f[i][j][7] = 0
                # ftemp[i_p][j_n][8] = f[i][j][8]#; f[i][j][8] = 0

    
    # ... and then computing macroscopic density and velocity for each lattice point        
    for j in range(1,sizeY_+1):
        for i in range(1,sizeX_+1):
            if solid[i][j] == 0:# or solid[i][j] == 1:
                rho[i][j] = 0
                ux[i][j] = 0
                uy[i][j] = 0
                for a in range(9):
                    rho[i][j] += ftemp[i][j][a]    
                    ux[i][j] += e_x[a]*ftemp[i][j][a]
                    uy[i][j] += e_y[a]*ftemp[i][j][a]
                ux[i][j] = ux[i][j]/rho[i][j] if rho[i][j] <> 0 else 0
                uy[i][j] = ux[i][j]/rho[i][j] if rho[i][j] <> 0 else 0

    #Computing equilibrium distribution function
    for j in range(1,sizeY_+ 1):
        for i in range(1,sizeX_+ 1):
            if solid[i][j] == 0:
                fct1 = w1*rho[i][j]
                fct2 = w2*rho[i][j]
                fct3 = w3*rho[i][j]

                uxeq[i][j] = ux[i][j]                                       #uxeq will incorporate external forces, if any
                uyeq[i][j] = uy[i][j] 

                uxsq = uxeq[i][j]*uxeq[i][j]
                uysq = uyeq[i][j]*uyeq[i][j]

                uxuy5 = uxeq[i][j] + uyeq[i][j]
                uxuy6 = -uxeq[i][j] + uyeq[i][j]
                uxuy7 = -uxeq[i][j] - uyeq[i][j]
                uxuy8 = uxeq[i][j] - uyeq[i][j]

                usq = uxsq + uysq

                feq[i][j][0] = fct1*(1.                              - f3*usq)
                feq[i][j][1] = fct2*(1. + f1*uxeq[i][j] + f2*uxsq    - f3*usq)
                feq[i][j][2] = fct2*(1. + f1*uyeq[i][j] + f2*uysq    - f3*usq)
                feq[i][j][3] = fct2*(1. - f1*uxeq[i][j] + f2*uxsq    - f3*usq)
                feq[i][j][4] = fct2*(1. - f1*uyeq[i][j] + f2*uysq    - f3*usq)
                feq[i][j][5] = fct3*(1. + f1*uxuy5 + f2*uxuy5*uxuy5  - f3*usq)
                feq[i][j][6] = fct3*(1. + f1*uxuy6 + f2*uxuy6*uxuy6  - f3*usq)
                feq[i][j][7] = fct3*(1. + f1*uxuy7 + f2*uxuy7*uxuy7  - f3*usq)
                feq[i][j][8] = fct3*(1. + f1*uxuy8 + f2*uxuy8*uxuy8  - f3*usq)

    #Collision step
    for j in range(1,sizeY_+1):
        for i in range(1,sizeX_+1):
            if solid[i][j] == 0:
                for a in range(9):
                    f[i][j][a] = ftemp[i][j][a] - (ftemp[i][j][a] - feq[i][j][a]) / tau
            # else:
            #     #if (is_interior(i,j)):
            #     # f[i+1][j][1] = ftemp[i][j][3]; f[i-1][j][3] = ftemp[i][j][1]
            #     # f[i][j+1][2] = ftemp[i][j][4]; f[i][j-1][4] = ftemp[i][j][2]
            #     # f[i+1][j+1][5] = ftemp[i][j][7]; f[i-1][j-1][7] = ftemp[i][j][5]
            #     # f[i-1][j+1][6] = ftemp[i][j][8]; f[i+1][j-1][8] = ftemp[i][j][6]
            #     dummy = f[i][j][1]; f[i][j][1] = f[i][j][3]; f[i][j][3] = dummy
            #     dummy = f[i][j][2]; f[i][j][2] = f[i][j][4]; f[i][j][4] = dummy
            #     dummy = f[i][j][5]; f[i][j][5] = f[i][j][7]; f[i][j][7] = dummy
            #     dummy = f[i][j][6]; f[i][j][6] = f[i][j][8]; f[i][j][8] = dummy

    densityM = zip(*rho)         #Transpose matrix rho
    print np.sum(densityM), np.sum(ux), np.sum(uy)
    ax = sns.heatmap(densityM, annot=False, vmin=0, vmax=3)
    plt.draw()
    plt.pause(0.1)
    plt.clf()
####################################################### OUTPUT ###########################################################################

        
