import numpy as np
import matplotlib.pyplot as plt

txt_files = ("p7.txt", "p7100.txt", "p71000.txt", "p710000.txt")
n_steps = (r"$n_{steps} = 10$", r"$n_{steps} = 10^2$", r"$n_{steps} = 10^3$", r"$n_{steps} = 10^4$")

# Plotting v(x)
for i in range(len(txt_files)):
    data = np.loadtxt(txt_files[i])

    x = data[:, 1] # x - vector
    v_x = data[:, 0] # v(x) values
    plt.plot(x, v_x, label = n_steps[i])

# Plotting u(x)
data = np.loadtxt('p2.txt')
x = data[:, 0] # x - vector
u_x = data[:, 1] # u(x) values
plt.plot(x, u_x, color = "black", linestyle = "dotted", label = r"$u(x)$")

# for all the plotparts
plt.title(r"$u(x)$ vs. $v(x)$")
plt.xlabel(r"$x\in[0, 1]$")
plt.ylabel(r"$v(x)$: aprroximated values to $u(x)$")
plt.ylim(-0.1, 1.85)
plt.legend()
plt.savefig("p7.pdf")
