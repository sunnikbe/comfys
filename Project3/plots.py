
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

%matplotlib qt

fig = plt.Figure()
ax = plt.axes(projection='3d')

r = np.loadtxt('RK4.txt')
T = 50
dt = 1e-3
n = int(T/dt)
t = np.linspace(0,T,n)

# 3D plot for M particles
M = 2
for i in range(M):
    ax.plot3D(r[:,0+3*i],r[:,1+3*i],r[:,2+3*i])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.figure()
for i in range(M):
    plt.plot(t,r[:,2+3*i])

