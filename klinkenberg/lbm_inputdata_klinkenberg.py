####################################################### PREPARATION ###########################################################################
from numpy import *
import scipy.io as sp 


# #importing grid
m = genfromtxt('../porestructure/klinkenberg_11b.dat', delimiter="\t")
m = m.transpose()

#name = "r45_u3"
name2 = "Re100_klinkenberg"

#What we have now:
# 1. channel_5.dat
# 2. channel_9.dat
# 3. channel_33.dat
# 4. cylinder_5.dat
# 5. cylinder_11.dat
# 6. cylinder_15.dat
# 7. cylinder_45.dat
# 8. cylinder2_45.dat

print m

sizeX_ = len(m) - 2         #length in x-direction
sizeY_ = len(m[0]) - 2        #length in y-direction

#The number of iteration
T = 300            #Total time used in the simulation
dt = 1             #time interval

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
Ft = zeros(T)

uxeq = zeros((sizeX_+2,sizeY_+2))                                       #uxeq will incorporate external forces, if any
uyeq = zeros((sizeX_+2,sizeY_+2)) 
uxsq = zeros((sizeX_+2,sizeY_+2))
uysq = zeros((sizeX_+2,sizeY_+2))
uxuy5 = zeros((sizeX_+2,sizeY_+2))
uxuy6 = zeros((sizeX_+2,sizeY_+2))
uxuy7 = zeros((sizeX_+2,sizeY_+2))
uxuy8 = zeros((sizeX_+2,sizeY_+2))
usq = zeros((sizeX_+2,sizeY_+2))
fct1 = zeros((sizeX_+2,sizeY_+2))
fct2 = zeros((sizeX_+2,sizeY_+2))                                       #uxeq will incorporate external forces, if any
fct3 = zeros((sizeX_+2,sizeY_+2))

ru = zeros(sizeY_+2)


#Constants used
Re_x = 100.
# ux0_left = 0.04
# ux0_right = ux0_left*1.5
rho0 = 0.8
rho_left = 1.
rho_right = 1.
#visc = ux0_left*r/Re_x             #kinematic viscosity
tau0 = 1.0 #3*visc + 0.5                
e_ = array([[0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0],[0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]])         
w = array([4.0/9.0, 1.0/9.0, 1.0/36.0])
c_eq = array([3., 9./2., 3./2.])
f_tol = 0.0001 

beta = 2.
gamma = log(beta)/2


####################################################### FUNCTIONS ###########################################################################
#To find out whether a m is in contact in fluid or not
def is_interior(i,j):       
    if m[i,j] == 1 and m[i-1,j] == 1 and m[i+1,j] == 1 and m[i,j-1] == 1 and m[i,j+1] == 1 and m[i-1,j-1] == 1 and m[i+1,j-1] == 1 and m[i-1,j+1] == 1 and m[i+1,j+1] == 1:
        return True

####################################################### INTIALIZATION ###########################################################################

####Initial Condition for density distribution, f
#Initialize density distribution f, ...

f_init = 1.
for j in range(1,sizeY_+ 1):
    for i in range(1, sizeX_+ 1):                 
        if m[i,j] == 0:
            for a in [0]:
                f[i,j,a] = (4./9.)*f_init 
            for a in [1,2,3,4]:
                f[i,j,a] = (1./9.)*f_init
            for a in [5,6,7,8]:
                f[i,j,a] = (1./36.)*f_init 

###Von Neumann Boundary condition
#Initializing flux boundary density distribution
for i in range(1,sizeX_+ 1):
    for j in range(1,sizeY_ + 1):
        if m[i,j] == 2:
            #West side
            #for a in [0,2,3,4,6,7]:
            f_init = rho_left
            for a in [0]:
                ftemp[i,j,a] = (4./9.)*f_init 
                f[i,j,a] = ftemp[i,j,a]
            for a in [2,3,4]:
                ftemp[i,j,a] = (1./9.)*f_init
                f[i,j,a] = ftemp[i,j,a]
            for a in [6,7]:
                ftemp[i,j,a] = (1./36.)*f_init
                f[i,j,a] = ftemp[i,j,a] 
        if m[i,j] == 4:
            #East side
            #for a in [0,1,2,4,5,8]:
            f_init = rho_right
            for a in [0]:
                ftemp[i,j,a] = (4./9.)*f_init
                f[i,j,a] = ftemp[i,j,a] 
            for a in [1,2,4]:
                ftemp[i,j,a] = (1./9.)*f_init
                f[i,j,a] = ftemp[i,j,a]
            for a in [5,8]:
                ftemp[i,j,a] = (1./36.)*f_init
                f[i,j,a] = ftemp[i,j,a] 

#print f[:,:,[0,2,3,4,6,7]]
###Reynold Number

#print "kinematic viscosity: ", visc
print "relaxation parameter, tau: ", tau0

