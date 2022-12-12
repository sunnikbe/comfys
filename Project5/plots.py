
import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

U = pa.cx_cube()
U.load('Output.bin')

V = np.loadtxt('Potential_map.txt')

N_xy = U.n_cols
n = U.n_slices
T = 0.002
h = 2.5e-5
x = np.linspace(0,1,N_xy)
y = np.linspace(0,1,N_xy)
t = np.linspace(0,T,n)

extent = [x[0],x[-1],y[0],y[-1]]

plt.figure()
plt.imshow(V,extent=extent,cmap=plt.get_cmap('viridis'))
plt.xlabel('x')
plt.ylabel('y')

for i in range(0,n,int((n-1)/2)):
    plt.figure()
    U_t = U[pa.single_slice,i]
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.abs(U_t.max()))
    plt.imshow(pa.abs(U_t),extent=extent, cmap=plt.get_cmap('viridis'), norm=norm)
    plt.text(0.95,0.95, f't = {t[i]*1e3:.1f}'+r'$\times10^3$', color='white',
                    horizontalalignment='right', verticalalignment='top', fontsize=12)
    plt.xlabel('x')
    plt.ylabel('y')

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
print(pa.cdot(P_y,P_y))

plt.figure()
plt.plot(y,np.abs(P_y))

plt.show()

# Probability over time, problem 7
prob = [pa.cdot(U[pa.single_slice, i], U[pa.single_slice, i]) for i in range(1, n)]
prob.insert(0, 1)
zeros = np.zeros(n)

# plt.title('Without slits') # ./main.exe 0 0
plt.title(r'With double slit, $v_0 = 1 \times 10^{10}$') # ./main.exe 2 1e10
plt.plot(t, 1 - np.abs(prob), label = r'$P = 1 - p_{ij}^n$')
plt.plot(t, zeros, label = r'$P = 1$')
plt.xlabel('time [s]')
plt.xticks(rotation = 25)
plt.ylabel(r'Deviation from normalized probability $(P = 1)$')
plt.legend()
plt.subplots_adjust(bottom = 0.18)
# plt.show()
# plt.savefig('probovertime.pdf') # ./main.exe 0 0
# plt.savefig('probovertimedslit.pdf') # ./main.exe 2 1e10
