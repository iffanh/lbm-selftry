####################################################### PREPARATION ###########################################################################
from numpy import *
import random
import scipy.io as sp 


# #importing grid
m = genfromtxt('../porestructure/empty_30.dat', delimiter="\t")
m = m.transpose()

#What we have now:
# 1. channel_5.dat
# 2. channel_9.dat
# 3. channel_33.dat
# 4. cylinder_5.dat
# 5. cylinder_11.dat
# 6. cylinder_15.dat
# 7. tank_10.dat

print m

sizeX_ = len(m) - 2         #length in x-direction
sizeY_ = len(m[0]) - 2        #length in y-direction

#The number of iteration
T = 20000            #Total time used in the simulation
dt = 1             #time interval

#m = [[0 for j in xrange(sizeY_ + 2)] for i in xrange(sizeX_ + 2)]                       #Presence of m or not, 1 means m

#Declaring variables

rho = zeros((sizeX_+2,sizeY_+2))                         #Density of the lattice point, 
ux = zeros((sizeX_+2,sizeY_+2))                          #Macroscopic velocity of the lattice point 
uy = zeros((sizeX_+2,sizeY_+2))
u = zeros((sizeX_+2,sizeY_+2))
uxeq = zeros((sizeX_+2,sizeY_+2))                        #Macroscopic velocity of the lattice point 
uyeq = zeros((sizeX_+2,sizeY_+2))
f = zeros((sizeX_+2,sizeY_+2,9))      #Density distribution of the a point f[x position][y position][index]
ftemp = zeros((sizeX_+2,sizeY_+2,9))
feq = zeros((sizeX_+2,sizeY_+2,9))
tau = zeros((sizeX_+2,sizeY_+2,9))
P = zeros((sizeX_+2,sizeY_+2))

uxeq = zeros((sizeX_+2,sizeY_+2))                                       #uxeq will incorporate external forces, if any
uyeq = zeros((sizeX_+2,sizeY_+2)) 
uxsq = zeros((sizeX_+2,sizeY_+2))
uysq = zeros((sizeX_+2,sizeY_+2))
uxuy5 = zeros((sizeX_+2,sizeY_+2))
uxuy6 = zeros((sizeX_+2,sizeY_+2))
uxuy7 = zeros((sizeX_+2,sizeY_+2))
uxuy8 = zeros((sizeX_+2,sizeY_+2))
usq = zeros((sizeX_+2,sizeY_+2))


#SCMP variables
psi = zeros((sizeX_+2,sizeY_+2))                        #Interaction potential
Fx = zeros((sizeX_+2,sizeY_+2))                        #Interaction potential
Fy = zeros((sizeX_+2,sizeY_+2))                        #Interaction potential


#Constants used
#Re_x = 1.
#ux0 = 1.
#uy0 = 0.2   
#visc = (ux0*1)/Re_x             #kinematic viscosity
tau0 = 2. #3*visc + 0.5                
e_ = array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],[0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])         
w = array([4.0/9.0, 1.0/9.0, 1.0/36.0])
c_eq = array([3., 9./2., 3./2.])
f_tol = 0.01

#SCMP variables
rho0 = 200.
psi0 = 4.
G = -120.
#Fx = 0.
#Fy = 0.

####################################################### FUNCTIONS ###########################################################################
#To find out whether a m is in contact in fluid or not
def is_interior(i,j):       
    if m[i,j] == 1 and m[i-1,j] == 1 and m[i+1,j] == 1 and m[i,j-1] == 1 and m[i,j+1] == 1 and m[i-1,j-1] == 1 and m[i+1,j-1] == 1 and m[i-1,j+1] == 1 and m[i+1,j+1] == 1:
        return True

####################################################### INTIALIZATION ###########################################################################

####Initial Condition for density distribution, f
#Initialize density distribution f, ...

f_init = 0.1
# for j in range(1,sizeY_+ 1):
#     for i in range(1, sizeX_+ 1):    
for j in range(sizeY_+ 2):
    for i in range(sizeX_+ 2):             
        if m[i,j] == 0:
            for a in range(9):
                f[i,j,a] = f_init + 0.01*random.uniform(0,1)

###Von Neumann Boundary condition
#Initializing flux boundary density distribution
# for i in range(1,sizeX_+ 1):
#     for j in range(1,sizeY_ + 1):
# for j in range(sizeY_+ 2):
#     for i in range(sizeX_+ 2): 
#         if m[i,j] == 2:
#             #West side
#             for a in [0,2,3,4,6,7]:
#                 f[i,j,a] = f_init

#print f[:,:,[0,2,3,4,6,7]]
###Reynold Number

#print "kinematic viscosity: ", visc
print "relaxation parameter, tau: ", tau0

