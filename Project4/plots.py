
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('Output.txt')
cycles = data[:,0]
epsilon = data[:,1]
m = data[:,2]
CV = data[:,3]
xi = data[:,4]

figure,axs = plt.subplots(2,2)
axs[0,0].plot(cycles,epsilon)
axs[0,1].plot(cycles,m)
axs[1,0].plot(cycles,CV)
axs[1,1].plot(cycles,xi)
plt.tight_layout()

plt.figure()
plt.hist(epsilon,bins=500)

plt.figure()
plt.hist(m,bins=500)

plt.show()

T = 1
N = 4
Z = 2*np.exp(8*T)+2*np.exp(-8*T) + 12
E = 1/Z*(-2*8*np.exp(8*T)+2*8*np.exp(-8*T))
E2 = 1/Z*(2*(-8)**2*np.exp(8*T)+2*8**2*np.exp(-8*T))
M = 1/Z*(2*4*np.exp(8*T))
M2 = 1/Z*(2*4**2*np.exp(8*T))

print(E,E2,M,M2)

eps = E/N
m = M/N

print(eps,m)

CV = 1/N/T**2*(E2 - E**2)
Xi = 1/N/T*(M2 - M**2)

print(CV,Xi)