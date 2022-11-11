
import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('Output.txt')
epsilon = data[:,0]
m = data[:,1]

print(np.mean(epsilon))
print(np.mean(m))

plt.figure()
plt.hist(epsilon,bins=500)

plt.figure()
plt.hist(m,bins=500)

plt.show()