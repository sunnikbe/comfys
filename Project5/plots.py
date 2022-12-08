
import numpy as np
import pyarma as pa
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import sys

slice = 0
if len(sys.argv) == 2:
    slice = int(sys.argv[1])

U = pa.cx_cube()
U.load('Output.bin')

V = np.loadtxt('Potential_map.txt')

N_xy = U.n_cols
N_t = U.n_slices

T = 0.002
h = 2.5e-5
n = int((T)/h)
x = np.linspace(0,1,N_xy)
y = np.linspace(0,1,N_xy)
t = np.linspace(0,T,n+1)

xmin,xmax = x[0],x[-1]
ymin,ymax = y[0],y[-1]

for i in range(0,n+1,int(n/2)):
    plt.figure()
    U_t = U[pa.single_slice,i]
    norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.abs(U_t.max()))
    plt.imshow(pa.abs(U_t),extent=[xmin,xmax,ymin,ymax], cmap=plt.get_cmap('viridis'), norm=norm)
    plt.text(0.95,0.95, f't = {t[i]*1e3:.2f} ms', color='white',
                    horizontalalignment='right', verticalalignment='top', fontsize=12)
    plt.xlabel('x')
    plt.ylabel('y')

fig = plt.figure()
ax = plt.gca()

U_0 = U[pa.single_slice,0]
norm = matplotlib.cm.colors.Normalize(vmin=0.0, vmax=np.abs(U_0.max()))

img = ax.imshow(pa.abs(U_0), extent=[xmin,xmax,ymin,ymax], cmap=plt.get_cmap('viridis'), norm=norm)

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
    time_txt.set_text(f't = {t[i]*1e3:.2f} ms')

    return img


anim = FuncAnimation(fig, animation, interval=250, frames=np.arange(0,U.n_slices,2))

anim.save("animation.mp4", writer="ffmpeg", bitrate=-1, fps=30)

plt.show()