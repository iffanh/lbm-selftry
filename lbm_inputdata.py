####################################################### PREPARATION ###########################################################################
from numpy import *
import scipy.io as sp 


# #importing grid
m = genfromtxt('porestructure/cylinder_11.dat', delimiter="\t")
m = m.transpose()

print m

sizeX_ = len(m) - 2         #length in x-direction
sizeY_ = len(m[0]) - 2        #length in y-direction

#The number of iteration
T = 1000            #Total time used in the simulation
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

#Constants used
tau0 = 1.2
e_ = array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],[0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])         
w = array([4.0/9.0, 1.0/9.0, 1.0/36.0])
c_eq = array([3., 9./2., 3./2.])
f_tol = 0.01 

####################################################### FUNCTIONS ###########################################################################
#To find out whether a m is in contact in fluid or not
def is_interior(i,j):       
    if m[i,j] == 1 and m[i-1,j] == 1 and m[i+1,j] == 1 and m[i,j-1] == 1 and m[i,j+1] == 1 and m[i-1,j-1] == 1 and m[i+1,j-1] == 1 and m[i-1,j+1] == 1 and m[i+1,j+1] == 1:
        return True

####################################################### INTIALIZATION ###########################################################################

####Initial Condition for density distribution, f
#Initialize density distribution f, ...

f_init = 0.1
for j in range(1,sizeY_+ 1):
    for i in range(1, sizeX_+ 1):                 
        if m[i,j] == 0:
            for a in range(9):
                f[i,j,a] = f_init 

###Von Neumann Boundary condition
#Initializing flux boundary density distribution
ux0 = 0.25
uy0 = 0.2

for i in range(1,sizeX_+ 1):
    for j in range(1,sizeY_ + 1):
        if m[i,j] == 2:
            #West side
            for a in [0,2,3,4,6,7]:
                f[i,j,a] = f_init

#print f[:,:,[0,2,3,4,6,7]]


