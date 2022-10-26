
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def convergence_rate_error(r_anal, r_anal_km1, r_num, r_num_km1, h_k, h_km1):
    delta_max_k = np.max(r_anal - r_num)
    delta_max_km1 = np.max(r_anal_km1 - r_num_km1)
    err = 1/3*np.log(delta_max_k/delta_max_km1)/np.log(h_k/h_km1)
    return err


def relative_error(r_anal,r_num):
    abs_err = np.linalg.norm(r_num - r_anal,axis=1)
    return (abs_err)/np.linalg.norm(r_anal,axis=1)


def solve_analytical(n_time_steps):
    t = np.linspace(0,50,n_time_steps)
    r_anal = np.zeros((3,n_time_steps))
    for i in range(n_time_steps):
        x = Ap*np.cos(omegap*t[i]) + Am*np.cos(omegam*t[i])
        y = -Ap*np.sin(omegap*t[i]) - Am*np.sin(omegam*t[i])
        z = z0*np.cos(omega_z*t[i])
        r_anal[:,i] = x,y,z


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
B0 = 96.5
V0d2 = 9.65
m = 40.078
omega_z2 = 2*q*V0d2/m
omega_z = np.sqrt(omega_z2)
omega0 = q*B0/m

omegap = (omega0 + np.sqrt(omega0**2-2*omega_z2))/2
omegam = (omega0 - np.sqrt(omega0**2-2*omega_z2))/2

Ap = (v0 + omegam*x0)/(omegam-omegap)
Am = -(v0 + omegap*x0)/(omegam-omegap)

x = Ap*np.cos(omegap*t) + Am*np.cos(omegam*t)
y = -Ap*np.sin(omegap*t) - Am*np.sin(omegam*t)
z = z0*np.cos(omega_z*t)

def r_analytically(N,T=50):
    t = np.linspace(0,T,N)
    x = Ap*np.cos(omegap*t) + Am*np.cos(omegam*t)
    y = -Ap*np.sin(omegap*t) - Am*np.sin(omegam*t)
    z = z0*np.cos(omega_z*t)
    r = np.zeros((N,3))
    for i in range(N):
        r[i,0] = x[i]
        r[i,1] = y[i]
        r[i,2] = z[i]
    return r,t


## Single particle test

r = np.loadtxt('RK4_pos_single_particle.txt')
r_euler = np.loadtxt('fEuler_pos.txt')

# plot z against time
plt.figure()
plt.plot(t,r[:,2])
plt.xlabel(r't ($\mu$s)')
plt.ylabel(r'z ($\mu$m)')
plt.savefig('plot0.pdf')

plt.figure()
N = np.array([4000,8000,16000,32000])

# relative error Runge-Kutta4
for N_i in N:
    h = T/N_i
    r_anal,t = r_analytically(N_i)
    r_num = np.loadtxt('RK4_error_calc_N_'+str(N_i)+'.txt')
    rel_err = relative_error(r_anal,r_num)
    plt.plot(t,rel_err,label=f'N = {N_i:.2e}')
plt.legend()
plt.xlabel(r't ($\mu$s)')
plt.ylabel(r'Relative error')
plt.savefig('plot1.pdf')

plt.figure()
# relative error forward Euler
for N_i in N:
    h = T/N_i
    r_anal,t = r_analytically(N_i)
    r_num = np.loadtxt('fEuler_error_calc_N_'+str(N_i)+'.txt')
    rel_err = relative_error(r_anal,r_num)
    plt.plot(t,rel_err,label=f'N = {N_i:.2e}')
plt.legend()
plt.xlabel(r't ($\mu$s)')
plt.ylabel(r'Relative error')
plt.savefig('plot2.pdf')

h = T/N
r_err_RK4 = 0
r_err_fEuler = 0
for i in range(1,4):
    r_anal_k, t = r_analytically(N[i])
    r_anal_km1, t = r_analytically(N[i-1])
    r_RK4_k = np.loadtxt('RK4_error_calc_N_'+str(N[i])+'.txt')
    r_RK4_km1 = np.loadtxt('RK4_error_calc_N_'+str(N[i-1])+'.txt')
    r_fEuler_k = np.loadtxt('fEuler_error_calc_N_'+str(N[i])+'.txt')
    r_fEuler_km1 = np.loadtxt('fEuler_error_calc_N_'+str(N[i-1])+'.txt')

    r_err_RK4 += convergence_rate_error(r_anal_k,r_anal_km1,r_RK4_k,r_RK4_km1,h[i],h[i-1])
    r_err_fEuler += convergence_rate_error(r_anal_k,r_anal_km1,r_fEuler_k,r_fEuler_km1,h[i],h[i-1])

print(f'r_err RK4 = {r_err_RK4:.2f}')
print(f'r_err fEuler = {r_err_fEuler:.2f}')

## Two particle tests

T = 50
dt = 1e-3
n = int(T/dt)
t = np.linspace(0,T,n)

# with Coulomb interaction
r = np.loadtxt('RK4_pos.txt')
v = np.loadtxt('RK4_vel.txt')

# without Coulomb interaction
r_no_int = np.loadtxt('RK4_pos_no_inter.txt')
v_no_int = np.loadtxt('RK4_vel_no_inter.txt')
M = 2 # number of particles

# plot x vs y
# with interactions
plt.figure()
for i in range(M):
    plt.plot(r[:,0+3*i],r[:,1+3*i])
plt.xlabel(r'x ($\mu$m)')
plt.ylabel(r'y ($\mu$m)')
plt.savefig('plot3.pdf')

# without interactions
plt.figure()
for i in range(M):
    plt.plot(r_no_int[:,0+3*i],r_no_int[:,1+3*i])
plt.xlabel(r'x ($\mu$m)')
plt.ylabel(r'y ($\mu$m)')
plt.savefig('plot4.pdf')

# phase plot x vs vx
# with interactions
plt.figure()
colors = ['red','blue']
for i in range(M):
    plt.plot(r[:,0+3*i],v[:,0+3*i],color=colors[i])
    plt.scatter(r[0,0+3*i],v[0,0+3*i],color=colors[i])
    plt.scatter(r[-1,0+3*i],v[-1,0+3*i],color=colors[i])
plt.xlabel(r'x ($\mu$m)')
plt.ylabel(r'v$_x$ ($\mu$m/$\mu$s)')
plt.axis('equal')
plt.savefig('plot5.pdf')

# without interactions
plt.figure()
for i in range(M):
    plt.plot(r_no_int[:,0+3*i],v_no_int[:,0+3*i],color=colors[i])
    plt.scatter(r_no_int[0,0+3*i],v_no_int[0,0+3*i],color=colors[i])
    plt.scatter(r_no_int[-1,0+3*i],v_no_int[-1,0+3*i],color=colors[i])
plt.xlabel(r'x ($\mu$m)')
plt.ylabel(r'v$_x$ ($\mu$m/$\mu$s)')
plt.axis('equal')
plt.savefig('plot6.pdf')

# phase plot z vs vz
# with interactions
plt.figure()
for i in range(M):
    plt.plot(r[:,2+3*i],v[:,2+3*i],color=colors[i])
    plt.scatter(r[0,2+3*i],v[0,2+3*i],color=colors[i])
    plt.scatter(r[-1,2+3*i],v[-1,2+3*i],color=colors[i])
plt.xlabel(r'z ($\mu$m)')
plt.ylabel(r'v$_z$ ($\mu$m/$\mu$s)')
plt.axis('equal')
plt.savefig('plot7.pdf')

# without interactions
plt.figure()
for i in range(M):
    plt.plot(r_no_int[:,2+3*i],v_no_int[:,2+3*i],color=colors[i])
    plt.scatter(r_no_int[0,2+3*i],v_no_int[0,2+3*i],color=colors[i])
    plt.scatter(r_no_int[-1,2+3*i],v_no_int[-1,2+3*i],color=colors[i])
plt.xlabel(r'z ($\mu$m)')
plt.ylabel(r'v$_z$ ($\mu$m/$\mu$s)')
plt.axis('equal')
plt.savefig('plot8.pdf')

plt.figure()
fig = plt.Figure()
ax = plt.axes(projection='3d')

# 3D plot for M particles with interactions
for i in range(M):
    ax.plot3D(r[:,0+3*i],r[:,1+3*i],r[:,2+3*i])

ax.set_xlabel(r'x ($\mu$m)')
ax.set_ylabel(r'y ($\mu$m)')
ax.set_zlabel(r'z ($\mu$m)')
plt.savefig('plot9.pdf')

plt.figure()
fig = plt.Figure()
ax = plt.axes(projection='3d')

# 3D plot for M particles without interactions
for i in range(M):
    ax.plot3D(r_no_int[:,0+3*i],r_no_int[:,1+3*i],r_no_int[:,2+3*i])

ax.set_xlabel(r'x ($\mu$m)')
ax.set_ylabel(r'y ($\mu$m)')
ax.set_zlabel(r'z ($\mu$m)')
plt.savefig('plot10.pdf')


# problem 9
plt.figure()
freq = np.array([0.1,0.4,0.7])
for i in range(len(freq)):
    data = np.loadtxt('part_in_trap_freq_'+str(freq[i])+'00000.txt')
    omega_v = data[:,0]
    part_in_trap = data[:,1]
    plt.plot(omega_v,part_in_trap/100,label=f'f = {freq[i]:.1f}')
plt.ylabel('Fraction of particles remaining')
plt.xlabel(r'$\omega_v$ (MHz)')
plt.legend()
plt.savefig('plot11.pdf')

plt.figure()
f = 0.4
data = np.loadtxt('C_i_resonans_trap.txt')
omega_v = data[:,0]
part_in_trap = data[:,1]
plt.plot(omega_v,part_in_trap)
plt.ylabel('Fraction of particles remaining')
plt.xlabel(r'$\omega_v$ (MHz)')
plt.savefig('plot12.pdf')

plt.figure()
f = 0.4
data = np.loadtxt('resonans_trap.txt')
omega_v = data[:,0]
part_in_trap = data[:,1]
plt.plot(omega_v,part_in_trap)
plt.ylabel('Fraction of particles remaining')
plt.xlabel(r'$\omega_v$ (MHz)')
plt.savefig('plot13.pdf')



plt.show()