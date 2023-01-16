
import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# terminal input for the number of slits used in the simulations
# only needed for filenames, won't matter for the plots
print('Number of slits: ')
number_of_slits = input()
print('End time T: ')
T = float(input())

# load the states from file generated by main.cpp
U = pa.cx_cube()
U.load('Output.bin')

# compute some variable used for the plots
N_xy = U.n_cols
n = U.n_slices
x = np.linspace(0,1,N_xy)
y = np.linspace(0,1,N_xy)
t = np.linspace(0,T,n)

# define the size of the slit simulation plots
extent = [x[0],x[-1],y[0],y[-1]]

# plot the absolute value of the state for the times t: 0.0, 0.001, 0.002
for i in range(0,n,int((n-1)/2)):
    plt.figure()
    U_t = U[pa.single_slice,i]
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.abs(U_t.max()))
    plt.imshow(pa.abs(U_t),extent=extent, cmap=plt.get_cmap('viridis'), norm=norm)
    plt.text(0.95,0.95, f't = {t[i]*1e3:.1f}'+r'$\times10^3$', color='white',
                    horizontalalignment='right', verticalalignment='top', fontsize=12)
    plt.xlabel('x')
    plt.ylabel('y')


# create a simulation of the wave packet for the double slit experiment
fig = plt.figure()
ax = plt.gca()

U_0 = U[pa.single_slice,0]
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.abs(U_0.max()))

img = ax.imshow(pa.abs(U_0), extent=extent, cmap=plt.get_cmap('viridis'), norm=norm)

ax.set_xlabel('x')
ax.set_ylabel('y')

cbar = fig.colorbar(img, ax=ax)
cbar.set_label('u(x,y,t)', fontsize=12)
cbar.ax.tick_params(labelsize=12)

time_txt = plt.text(0.95,0.95, f't = {t[0]*1e3:.2f} ms', color='white',
                    horizontalalignment='right', verticalalignment='top', fontsize=12)

def animation(i):
    U_i = U[pa.single_slice,i]
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.abs(U_i.max()))
    img.set_norm(norm)

    img.set_data(pa.abs(U_i))
    time_txt.set_text(f't = {t[i]*1e3:.2f}'+r'$\times10^3$')

    return img


anim = FuncAnimation(fig, animation, interval=250, frames=np.arange(0,U.n_slices,2))

anim.save("animation.mp4", writer="ffmpeg", bitrate=-1, fps=30)


# probability in y direction of the wave packet along x = 0.8 at time t = 0
# for single-, double- and triple-slits problem
U_t = pa.normalise(U[pa.single_slice,n-1],2,0)
xi = int(0.8*N_xy)
P_y = U_t[:,xi]

plt.figure()
plt.plot(y,np.abs(P_y))
plt.xlabel('y')
plt.ylabel('P(y|x=0.8;t=0.002')
plt.savefig(f'probplot_{number_of_slits}_slits.pdf')

# Probability over time, problem 7
prob = [pa.cdot(U[pa.single_slice, i], U[pa.single_slice, i]) for i in range(1, n)]
prob.insert(0, 1)
zeros = np.zeros(n)

plt.figure()
plt.title('Without slits') # ./main.exe 0 0 0.05 0.008
#plt.title(r'With double slit, $v_0 = 1 \times 10^{10}$') # ./main.exe 2 1e10 0.1 0.008
plt.plot(t, 1 - np.abs(prob), label = r'$P = 1 - p_{ij}^n$')
plt.plot(t, zeros, label = r'$P = 1$')
plt.xlabel('t')
plt.xticks(rotation = 25)
plt.ylabel(r'Deviation from normalized probability $(P = 1)$')
plt.legend()
plt.subplots_adjust(bottom = 0.18)
plt.savefig('probovertime.pdf') # ./main.exe 0 0 0.05 0.008
#plt.savefig('probovertimedslit.pdf') # ./main.exe 2 1e10 0.1 0.008

# Problem 8 Re(u_ij) Im(u_ij)
for i in range(0,n,int((n-1)/2)):
    plt.figure()
    U_t = U[pa.single_slice, i]
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.abs(U_t.max()))
    plt.imshow(pa.imag(U_t),extent=extent, cmap=plt.get_cmap('viridis'), norm=norm)
    plt.text(0.95,0.95, f't = {t[i]*1e3:.1f}'+r'$\times10^3$', color='white',
                    horizontalalignment='right', verticalalignment='top', fontsize=12)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'Im($u_{ij}$)')

# Re(u_ij)
for i in range(0,n,int((n-1)/2)):
    plt.figure()
    U_t = U[pa.single_slice, i]
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.abs(U_t.max()))
    plt.imshow(pa.real(U_t),extent=extent, cmap=plt.get_cmap('viridis'), norm=norm)
    plt.text(0.95,0.95, f't = {t[i]*1e3:.1f}'+r'$\times10^3$', color='white',
                    horizontalalignment='right', verticalalignment='top', fontsize=12)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title(r'Re($u_{ij}$)')

plt.show()
