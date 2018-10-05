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
tau = 1.5
e_x = [0.0, 1.0, 0.0, -1.0, 0.0, 1.0, -1.0, -1.0, 1.0]          
e_y = [0.0, 0.0, 1.0, 0.0, -1.0, 1.0, 1.0, -1.0, -1.0]
w = [4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0]
w1 = 4./9.
w2 = 1./9.
w3 = 1./36.
f1 = 3.
f2 = 9./2.
f3 = 3./2.
f_tol = 0. 

####################################################### FUNCTIONS ###########################################################################
#To find out whether a solid is in contact in fluid or not
def is_interior(i,j):       
    if solid[i][j] == 1 and solid[i-1][j] == 1 and solid[i+1][j] == 1 and solid[i][j-1] == 1 and solid[i][j+1] == 1 and solid[i-1][j-1] == 1 and solid[i+1][j-1] == 1 and solid[i-1][j+1] == 1 and solid[i+1][j+1] == 1:
        return True

####################################################### INTIALIZATION ###########################################################################

####Sacrificed grid 
for i in range(0,sizeX_+ 2):
    for j in range(0,sizeY_ + 2):
        solid[i][0] = 0
        solid[i][sizeY_+1] = 0
        solid[0][j] = 0
        solid[sizeX_+1][j] = 0

####Solid presence in Grid 
#The outer boundary
#0 is a normal no solid boundary
#1 is a bounce back boundary
#2 is a von neumann boundary
for i in range(1,sizeX_+ 1):
    for j in range(1,sizeY_ + 1):
        solid[i][1] = 1
        solid[i][sizeY_] = 1
        solid[1][j] = 2
        solid[sizeX_][j] = 0

solid[1][1] = 1
solid[1][sizeY_] = 1
solid[sizeX_][1] = 1
solid[sizeX_][sizeY_] = 1

#The cylinder
solid[7][10] = 1; solid[8][10] = 1; solid[9][10] = 1
solid[6][9] = 1; solid[7][9] = 1; solid[8][9] = 1; solid[9][9] = 1; solid[10][9] = 1
solid[6][8] = 1; solid[7][8] = 1; solid[8][8] = 1; solid[9][8] = 1; solid[10][8] = 1
solid[6][7] = 1; solid[7][7] = 1; solid[8][7] = 1; solid[9][7] = 1; solid[10][7] = 1
solid[7][6] = 1; solid[8][6] = 1; solid[9][6] = 1

####Initial Condition for density distribution, f
#Initialize density distribution f, ...

f_init = 0.1
for j in range(1,sizeY_+ 1):
    for i in range(1, sizeX_+ 1):                 
        if solid[i][j] == 0:
            for a in range(9):
                f[i][j][a] = f_init 

###Von Neumann Boundary condition
#Initializing flux boundary density distribution
ux0 = 0.8
for i in range(1,sizeX_+ 1):
    for j in range(1,sizeY_ + 1):
        if solid[i][j] == 2:
            #West side
            for a in [0,2,3,4,6,7]:
                f[i][j][a] = f_init

#f[10][10][5] = 10.