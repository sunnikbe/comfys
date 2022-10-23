
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#%matplotlib qt

## variables used in Penning trap

T = 50
dt = 1e-3
n = int(T/dt)
t = np.linspace(0,T,n)

## analytical solutions

x0 = 20
v0 = 25
z0 = 20
q = 1
V0d2 = 9.65
m = 40.078
omega_z = np.sqrt(2*q*V0d2/m)
z = z0*np.cos(omega_z*t)


## Single particle test

r = np.loadtxt('RK4_pos_single_particle.txt')
v = np.loadtxt('RK4_vel_single_particle.txt')
r_euler = np.loadtxt('fEuler_pos.txt')

# plot z against time
plt.figure()
plt.plot(t,r[:,2])
plt.plot(t,z)
plt.plot(t,r_euler[:,2])

plt.figure()
plt.plot(t,r[:,0])

plt.figure()
plt.plot(t,r[:,1])

plt.figure()
fig = plt.Figure()
ax = plt.axes(projection='3d')

# 3D plot for M particles
ax.plot3D(r[:,0],r[:,1],r[:,2])
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.figure()
fig = plt.Figure()
ax = plt.axes(projection='3d')

# 3D plot for M particles
ax.plot3D(r_euler[:,0],r_euler[:,1],r_euler[:,2])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

## Two particle tests

r = np.loadtxt('RK4_pos.txt')
v = np.loadtxt('RK4_vel.txt')
M = 2

plt.figure()
for i in range(M):
    plt.plot(r[:,0+3*i],r[:,1+3*i])

plt.figure()
for i in range(M):
    plt.plot(r[:,0+3*i],v[:,0+3*i])
plt.xlabel('x')
plt.ylabel('v_x')
plt.axis('equal')

plt.figure()
for i in range(M):
    plt.plot(r[:,2+3*i],v[:,2+3*i])
plt.xlabel('z')
plt.ylabel('v_z')
plt.axis('equal')

plt.figure()
fig = plt.Figure()
ax = plt.axes(projection='3d')

# 3D plot for M particles
for i in range(M):
    ax.plot3D(r[:,0+3*i],r[:,1+3*i],r[:,2+3*i])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')



plt.show()