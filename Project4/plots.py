
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

def read_from_file(L,T,filename,n_threads=1):
    N = L**2
    if n_threads == 1:
        E = np.cumsum(np.loadtxt('E_mat_'+filename+'.txt'))
        EE = np.cumsum(np.loadtxt('EE_mat_'+filename+'.txt'))
        M = np.cumsum(np.loadtxt('M_mat_'+filename+'.txt'))
        MM = np.cumsum(np.loadtxt('MM_mat_'+filename+'.txt'))
    else:
        E = np.cumsum(np.sum(np.loadtxt('E_mat_'+filename+'.txt'),axis=1)/n_threads)
        EE = np.cumsum(np.sum(np.loadtxt('EE_mat_'+filename+'.txt'),axis=1)/n_threads)
        M = np.cumsum(np.sum(np.loadtxt('M_mat_'+filename+'.txt'),axis=1)/n_threads)
        MM = np.cumsum(np.sum(np.loadtxt('MM_mat_'+filename+'.txt'),axis=1)/n_threads)
    n = E.shape[0]
    cycles = np.arange(1,n+1)

    eps = np.zeros(n)
    m = np.zeros(n)
    CV = np.zeros(n)
    xi = np.zeros(n)
    for i in range(1,n+1):
        eps[i-1] = E[i-1]/N/i
        m[i-1] = M[i-1]/N/i
        CV[i-1] = 1/N/T**2*(EE[i-1]/i - (E[i-1]/i)**2)
        xi[i-1] = 1/N/T*(MM[i-1]/i - (M[i-1]/i)**2)
    return eps,m,CV,xi,cycles


## Numerical solutions

# 2x2 lattice case with T = 1.0
eps,m,CV,xi,cycles = read_from_file(2,1.0,"set1")
figure,axs = plt.subplots(2,2)
axs[0,0].plot(cycles,eps)
axs[0,1].plot(cycles,m)
axs[1,0].plot(cycles,CV)
axs[1,1].plot(cycles,xi)
axs[0,0].set_xlabel('cycles')
axs[0,0].set_ylabel(r'$\epsilon$ (J)')
axs[0,1].set_xlabel('cycles')
axs[0,1].set_ylabel(r'|m|')
axs[1,0].set_xlabel('cycles')
axs[1,0].set_ylabel(r'$C_V$ ($k_B$)')
axs[1,1].set_xlabel('cycles')
axs[1,1].set_ylabel(r'$\chi$')
plt.tight_layout()
plt.savefig('2x2_exp_values.pdf')

L = 20
T = [1.0,1.0,2.4,2.4]
filename = ['set1_0_ordered','set1_0_unordered','set2_4_ordered','set2_4_unordered']
label = ['ordered','unordered']
colors = ['red','blue','green','purple']
lines = ['--',':']
figure,axs = plt.subplots(2,1)
for i in range(4):
    eps,m,CV,xi,cycles = read_from_file(L,T[i],filename[i])
    axs[0].plot(cycles,eps,linestyle=lines[i%2],color=colors[i],label=f'T = {T[i]:.1f} '+r'J/k$_B$'+f', {label[i%2]}')
    axs[1].plot(cycles,m,linestyle=lines[i%2],color=colors[i],label=f'T = {T[i]:.1f} '+r'J/k$_B$'+f', {label[i%2]}')
axs[0].legend(loc='best',bbox_to_anchor=(1.02,1))
axs[1].legend(loc='best',bbox_to_anchor=(1.02,1))
axs[0].set_xlabel('cycles')
axs[0].set_ylabel(r'$\epsilon$ (J)')
axs[1].set_xlabel('cycles')
axs[1].set_ylabel(r'|m|')
plt.tight_layout()
plt.savefig('20x20_exp_values.pdf')

x_min = [-1.9978,-1.3]
x_max = [-1.99675,-1.22]
figure,axs = plt.subplots(1,2)
for i in range(2):
    eps,m,CV,xi,cycles = read_from_file(L,T[i],filename[2*i])
    sns.histplot(eps, kde=True, bins=500, ax=axs[i], stat='probability')
    axs[i].set_xlim(x_min[i],x_max[i])
plt.savefig('20x20_exp_values_histogram.pdf')

L = [40,60,80,100]
figure,axs = plt.subplots(2,2)
for Ls in L:
    data = np.loadtxt('crit_temp_L_'+str(Ls)+'.txt')
    T = data[:,0]
    eps = data[:,1]
    m = data[:,2]
    CV = data[:,3]
    xi = data[:,4]

    axs[0,0].plot(T,eps)
    axs[0,1].plot(T,m)
    axs[1,0].plot(T,CV)
    axs[1,1].plot(T,xi)

axs[0,0].set_xlabel(r'T (J/$k_B$)')
axs[0,0].set_ylabel(r'$\epsilon$ (J)')
axs[0,1].set_xlabel(r'T (J/$k_B$)')
axs[0,1].set_ylabel(r'|m|')
axs[1,0].set_xlabel(r'T (J/$k_B$)')
axs[1,0].set_ylabel(r'$C_V$ ($k_B$)')
axs[1,1].set_xlabel(r'T (J/$k_B$)')
axs[1,1].set_ylabel(r'$\chi$')
plt.tight_layout()
plt.savefig('phase_transitions.pdf')

"""T_c_L = [2.3,]
lin_reg = stats.linregress(1/L,T_c_L)
alpha = lin_reg.slope
T_c = lin_reg.intercept
plt.plot(1/L,alpha*T_c_L+T_c)"""

## Analytical computations

T = 1.0
N = 4
Z = 2*np.exp(8/T)+2*np.exp(-8/T) + 12
E = 1/Z*(-2*8*np.exp(8/T)+2*8*np.exp(-8/T)) 
E2 = 1/Z*(2*(-8/T)**2*np.exp(8/T)+2*8**2*np.exp(-8/T))
M = 1/Z*(2*4*np.exp(8/T)+8*2)
M2 = 1/Z*(2*4**2*np.exp(8/T)+8*2**2)

print(E,E2,M,M2)

eps = E/N
m = M/N

print(eps,m)

CV = (E2 - E**2)/N/T/T
Xi = (M2 - M**2)/N/T

print(CV,Xi)



plt.show()