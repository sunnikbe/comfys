
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

# read a file and compute the energy and magnetization per spin,
# the heat capacity and susceptiblity
def read_from_file(L,T,filename):
    N = L**2
    E = np.cumsum(np.loadtxt('E_mat_'+filename+'.txt'))/N
    EE = np.cumsum(np.loadtxt('EE_mat_'+filename+'.txt'))/N
    M = np.cumsum(np.loadtxt('M_mat_'+filename+'.txt'))/N
    MM = np.cumsum(np.loadtxt('MM_mat_'+filename+'.txt'))/N
    n = E.shape[0]
    cycles = np.arange(1,n+1)

    eps = np.zeros(n)
    m = np.zeros(n)
    CV = np.zeros(n)
    chi = np.zeros(n)
    for i in range(1,n+1):
        eps[i-1] = E[i-1]/N/i
        m[i-1] = M[i-1]/N/i
        CV[i-1] = 1/N/T**2*(EE[i-1]/i - (E[i-1]/i)**2)
        chi[i-1] = 1/N/T*(MM[i-1]/i - (M[i-1]/i)**2)
    return eps,m,CV,chi,cycles


## Analytical computations

# solve the 2x2 size lattive Ising model analytically for a given temperature
T = 1.0
N = 4
Z = 2*np.exp(8/T)+2*np.exp(-8/T) + 12
E = 1/Z*(-2*8*np.exp(8/T)+2*8*np.exp(-8/T)) 
E2 = 1/Z*(2*(-8/T)**2*np.exp(8/T)+2*8**2*np.exp(-8/T))
M = 1/Z*(2*4*np.exp(8/T)+8*2)
M2 = 1/Z*(2*4**2*np.exp(8/T)+8*2**2)

eps_anal = E/N
m_anal = M/N

print(f'epsilon = {eps_anal:.3f} J')
print(f'm = {m_anal:.3f}')

CV_anal = (E2 - E**2)/N/T/T
chi_anal = (M2 - M**2)/N/T

print(f'C_V = {CV_anal:.3f} k_B')
print(f'chi = {chi_anal:.3f}')

## Numerical solutions

plt.rcParams['font.size'] = '12'

# 2x2 lattice case with T = 1.0
# plots the energy and magnetization per spin, heat capacity and susceptebility
# together with the analytical solution
eps,m,CV,chi,cycles = read_from_file(2,1.0,"set1")
figure,axs = plt.subplots(2,2)
axs[0,0].plot(cycles,eps,label='Numerical')
axs[0,0].plot(np.linspace(eps_anal,eps_anal,len(eps)),'--',label='Analytical')
axs[0,1].plot(cycles,m,label='Numerical')
axs[0,1].plot(np.linspace(m_anal,m_anal,len(m)),'--',label='Analytical')
axs[1,0].plot(cycles,CV,label='Numerical')
axs[1,0].plot(np.linspace(CV_anal,CV_anal,len(CV)),'--',label='Analytical')
axs[1,1].plot(cycles,chi,label='Numerical')
axs[1,1].plot(np.linspace(chi_anal,chi_anal,len(chi)),'--',label='Analytical')
axs[0,0].legend(fontsize=9)
axs[1,0].legend(fontsize=9)
axs[0,1].legend(fontsize=9)
axs[1,1].legend(fontsize=9)
axs[0,0].set_xlabel('cycles')
axs[0,0].set_ylabel(r'<$\epsilon$> (J)')
axs[0,1].set_xlabel('cycles')
axs[0,1].set_ylabel(r'<|m|>')
axs[1,0].set_xlabel('cycles')
axs[1,0].set_ylabel(r'$C_V$ ($k_B$)')
axs[1,1].set_xlabel('cycles')
axs[1,1].set_ylabel(r'$\chi$')
plt.tight_layout()
plt.savefig('2x2_exp_values.pdf')

# plot the 20x20 size lattice for T = 1.0 and T = 2.4 with both unordered and ordered
# starting spin configurations
L = 20
T = [1.0,1.0,2.4,2.4]
filename = ['set1_0_ordered','set1_0_unordered','set2_4_ordered','set2_4_unordered']
label = ['ordered','unordered']
colors = ['red','blue','green','purple']
lines = ['--',':']
figure,axs = plt.subplots(2,1)
for i in range(4):
    eps,m,CV,chi,cycles = read_from_file(L,T[i],filename[i])
    axs[0].plot(cycles,eps,linestyle=lines[i%2],color=colors[i],label=f'T = {T[i]:.1f} '+r'J/k$_B$'+f', {label[i%2]}')
    axs[1].plot(cycles,m,linestyle=lines[i%2],color=colors[i],label=f'T = {T[i]:.1f} '+r'J/k$_B$'+f', {label[i%2]}')
axs[0].legend(loc='best',bbox_to_anchor=(1.02,1))
axs[1].legend(loc='best',bbox_to_anchor=(1.02,1))
axs[0].set_xlabel('cycles')
axs[0].set_ylabel(r'<$\epsilon$> (J)')
axs[1].set_xlabel('cycles')
axs[1].set_ylabel(r'<|m|>')
plt.tight_layout()
plt.savefig('20x20_exp_values.pdf')

# plot a histogram of the sampled energy per spin for T = 1,0 and T = 2.4 
# with one million MCMC cycles
filename = ['hist_1_0','hist_2_4']
figure,axs = plt.subplots(1,2)
N = 400
for i in range(2):
    E = np.loadtxt('E_mat_'+filename[i]+'.txt')/N
    eps = E/N
    sns.histplot(eps, bins=500, ax=axs[i], stat='probability')
axs[0].set_xlabel(r'<$\epsilon$> (J)')
axs[1].set_xlabel(r'<$\epsilon$> (J)')
plt.tight_layout()
plt.savefig('20x20_exp_values_histogram.pdf')

# plot the energy and magnetization per spin, heat capacity and susceptibility
# for temperatures between 2.1 and 2.4
L = [40,60,80,100]
figure,axs = plt.subplots(2,2)
for Ls in L:
    data = np.loadtxt('crit_temp_L_'+str(Ls)+'.txt')
    T = data[:,0]
    eps = data[:,1]
    m = data[:,2]
    CV = data[:,3]
    chi = data[:,4]

    axs[0,0].plot(T,eps,label=f'L = {Ls}')
    axs[0,1].plot(T,m,label=f'L = {Ls}')
    axs[1,0].plot(T,CV,label=f'L = {Ls}')
    axs[1,1].plot(T,chi,label=f'L = {Ls}')

axs[0,0].set_xlabel(r'T (J/$k_B$)')
axs[0,0].set_ylabel(r'<$\epsilon$> (J)')
axs[0,1].set_xlabel(r'T (J/$k_B$)')
axs[0,1].set_ylabel(r'<|m|>')
axs[1,0].set_xlabel(r'T (J/$k_B$)')
axs[1,0].set_ylabel(r'$C_V$ ($k_B$)')
axs[1,1].set_xlabel(r'T (J/$k_B$)')
axs[1,1].set_ylabel(r'$\chi$')
axs[0,0].legend(fontsize=8)
axs[1,0].legend(fontsize=8)
axs[0,1].legend(fontsize=8)
axs[1,1].legend(fontsize=8)
plt.tight_layout()
plt.savefig('phase_transitions.pdf')

# create a linear regression for the critical temperatures found from the phase transition
# as a function of the inverse lattice size
plt.figure()
T_c_L = np.array([2.303,2.304,2.293,2.295])
L = np.array([40,60,80,100])
lin_reg = stats.linregress(1.0/L,T_c_L)
alpha = lin_reg.slope
T_c = lin_reg.intercept
print(f'T_c = {T_c:.3f} J/k_B')
plt.plot(1/L,alpha*1/L+T_c)
plt.xlabel('1/L')
plt.ylabel(r'T$_c$ (J/k_B)')
plt.savefig('linear_regression.pdf')


plt.show()