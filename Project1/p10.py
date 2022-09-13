import numpy as np

# General algorithm
g_data = np.loadtxt("p10.txt")
t_10 = g_data[:5, 1]
t_102 = g_data[5:10, 1]
t_103 = g_data[10:15, 1]
t_104 = g_data[15:20, 1]
t_105 = g_data[20:25, 1]
t_106 = g_data[25:30, 1]

# Special algorithm
s_data = np.loadtxt("p10s.txt")
t_10_s = s_data[:5, 1]
t_102_s = s_data[5:10, 1]
t_103_s = s_data[10:15, 1]
t_104_s = s_data[15:20, 1]
t_105_s = s_data[20:25, 1]
t_106_s = s_data[25:30, 1]

#Finding means and stderr = sqrt(std^2)
t_g = [t_10, t_102, t_103, t_104, t_105, t_106]
t_s = [t_10_s, t_102_s, t_103_s, t_104_s, t_105_s, t_106_s]

meanvals_t_g = []
meanvals_t_s = []

stderr_t_g = []
stderr_t_s = []

print('General alg')

for i in range(len(t_g)):
    #general
    meanvals_t_g.append(np.mean(t_g[i]))
    stderr_t_g.append(np.sqrt(np.std(t_g[i])**2))

    #special
    meanvals_t_s.append(np.mean(t_s[i]))
    stderr_t_s.append(np.sqrt(np.std(t_s[i])**2))

    print(f'{meanvals_t_g[i]:.2E} +- {stderr_t_g[i]:.2E}')

print('Special alg')
for j in range(len(t_s)):
    print(f'{meanvals_t_s[j]:.2E} +- {stderr_t_s[j]:.2E}')


"""
Terminal:

Sunnivas-MacBook-Pro:Project1 sunnivakbergan$ python3 p10.py
General alg
4.20E-06 +- 4.00E-07
8.60E-06 +- 1.36E-06
5.20E-05 +- 9.51E-06
4.86E-04 +- 7.27E-05
5.04E-03 +- 7.45E-04
4.86E-02 +- 7.38E-03
Special alg
4.60E-06 +- 1.20E-06
6.20E-06 +- 9.80E-07
3.20E-05 +- 2.53E-06
2.79E-04 +- 8.22E-06
3.31E-03 +- 7.79E-04
3.18E-02 +- 6.89E-03
"""
