import matplotlib.pyplot as plt
import numpy as np 
import os

#Plotting dragforce in the cylinder
t = range(1000)
a = np.load(os.path.join("r5_u3", "dragForce_r5_u3.npy"))
b = np.load(os.path.join("r11_u3", "dragForce_r11_u3.npy"))
c = np.load(os.path.join("r15_u3", "dragForce_r15_u3.npy"))
d = np.load(os.path.join("r45_u3", "dragForce_r45_u3.npy"))

plt.figure(1)
plt.plot(t, a, label='r1')
plt.plot(t, b, label='r2')
plt.plot(t, c, label='r3')
plt.plot(t, d, label='r4')
plt.legend(loc='upper right')
plt.show()