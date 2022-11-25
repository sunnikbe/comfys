
import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt('Potential_map.txt')
N = data.shape[0]
x = np.linspace(0,1,N)
y = np.linspace(0,1,N)

plt.contourf(x,y,data)
plt.xlabel('x')
plt.ylabel('y')
plt.colorbar()



plt.show()