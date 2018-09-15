#This is an attempt to make a program to simulate fluid flow using LBM
#This is for single-phase flow

import math
####################################################### PREPARATION ###########################################################################
#The size of grid

sizeX_ = 5         #length in x-direction
sizeY_ = 4         #length in y-direction

solid = [[0 for j in xrange(sizeY_)] for i in xrange(sizeX_)]                       #Presence of solid or not, 1 means solid

#Declaring variables

rho = [[0 for j in xrange(sizeY_)] for i in xrange(sizeX_)]                         #Density of the lattice point, 
ux = [[0 for j in xrange(sizeY_)] for i in xrange(sizeX_)]                               #Macroscopic velocity of the lattice point 
uy = [[0 for j in xrange(sizeY_)] for i in xrange(sizeX_)]
uxeq = [[0 for j in xrange(sizeY_)] for i in xrange(sizeX_)]                               #Macroscopic velocity of the lattice point 
uyeq = [[0 for j in xrange(sizeY_)] for i in xrange(sizeX_)]
f = [[[0 for k in xrange(9)] for j in xrange(sizeY_)] for i in xrange(sizeX_)]      #Density distribution of the a point f[x position][y position][index]
ftemp = [[[0 for k in xrange(9)] for j in xrange(sizeY_)] for i in xrange(sizeX_)]
feq = [[[0 for k in xrange(9)] for j in xrange(sizeY_)] for i in xrange(sizeX_)]

#Constants used
e_x = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0]          
e_y = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]
w = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]
w1 = 4./9.
w2 = 1./9.
w3 = 1./36.
f1 = 3.
f2 = 9./2.
f3 = 3./2. 

####################################################### INTIALIZATION ###########################################################################

####Initial Condition

#Initialize density distribution f, 
for j in range(sizeY_):
    for i in range(sizeX_):
        for a in range(9):
            f[i][j][a] = 0.8 if i < (sizeX_/2) else 0.1 
            
# and then computing macroscopic density and velocity for each lattice point        
for j in range(sizeY_):
    for i in range(sizeX_):
        for a in range(9):
            if solid[i][j] == 0:
                rho[i][j] += f[i][j][a]
                ux[i][j] += e_x[a]*f[i][j][a]
                uy[i][j] += e_y[a]*f[i][j][a]
            ux[i][j] /= rho[i][j]
            uy[i][j] /= rho[i][j]
            

#Streaming step
for j in range(0,sizeY_):
    j_n = (j-1) if j > 0 else (sizeY_-1)
    j_p = (j+1) if j < (sizeY_ - 1) else 0 

    for i in range(0,sizeX_):
        i_n = (i-1) if i > 0 else (sizeX_-1)
        i_p = (i+1) if i < (sizeX_ - 1) else 0 

        if solid[i][j] == 0:
            ftemp[i][j][0] = f[i][j][0]
            ftemp[i_p][j][1] = f[i][j][1]
            ftemp[i][j_p][2] = f[i][j][2]
            ftemp[i_n][j][3] = f[i][j][3]
            ftemp[i][j_n][4] = f[i][j][4]
            ftemp[i_p][j_p][5] = f[i][j][5]
            ftemp[i_n][j_p][6] = f[i][j][6]
            ftemp[i_n][j_n][7] = f[i][j][7]
            ftemp[i_p][j_n][8] = f[i][j][8]


#Computing equilibrium distribution function
for j in range(0,sizeY_):
    for i in range(0,sizeX_):
        if solid[i][j] == 0:
            fct1 = w1*rho[i][j]
            fct2 = w2*rho[i][j]
            fct3 = w3*rho[i][j]

            uxeq[i][j] = ux[i][j]                                       #exeq will incorporate external forces, if any
            uyeq[i][j] = uy[i][j] 

            uxsq = uxeq[i][j]*uxeq[i][j]
            uysq = uyeq[i][j]*uyeq[i][j]

            uxuy5 = uxeq[i][j] + uyeq[i][j]
            uxuy6 = -uxeq[i][j] + uyeq[i][j]
            uxuy7 = -uxeq[i][j] + uyeq[i][j]
            uxuy8 = uxeq[i][j] - uyeq[i][j]

            usq = uxsq + uysq

            feq[i][j][0] = fct1*(1.                              - f3*usq)
            feq[i][j][1] = fct2*(1. + f1*uxeq[i][j] + f2*uxsq    - f3*usq)
            feq[i][j][2] = fct2*(1. + f1*uxeq[i][j] + f2*uxsq    - f3*usq)
            feq[i][j][3] = fct2*(1. - f1*uxeq[i][j] + f2*uxsq    - f3*usq)
            feq[i][j][4] = fct2*(1. - f1*uxeq[i][j] + f2*uxsq    - f3*usq)
            feq[i][j][5] = fct3*(1. + f1*uxuy5 + f2*uxuy5*uxuy5  - f3*usq)
            feq[i][j][6] = fct3*(1. + f1*uxuy6 + f2*uxuy6*uxuy6  - f3*usq)
            feq[i][j][7] = fct3*(1. + f1*uxuy7 + f2*uxuy7*uxuy7  - f3*usq)
            feq[i][j][8] = fct3*(1. + f1*uxuy8 + f2*uxuy8*uxuy8  - f3*usq)


print feq
print ftemp

#Feq value is much higher than f value, try to fix it!            
