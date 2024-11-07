"""
Relative permeability, as a function of pressure
Van Genuchten specification, Ref [VG80]
cf PF-UM

Created on Aug 28 2024, 3:41pm
author : SP.
"""

import matplotlib.pyplot as plt
import numpy as np

plt.style.use('/home/patras/PF-Valmalenco/scripts/config.mplstyle')
#path_fig = '/mnt/c/Users/User/Documents/POLIMI/0_TESI/8-Figures/'   #local
path_fig = '/mnt/c/Users/Sophie/Documents/4-Figures/'   #distant

alpha_val = [1.0,1.0,1.8]
n_val = [1.4,3.0,1.3]
nb = len(alpha_val)

# VanGenuchten
def kr(p,alpha,n):
    m = 1-1/n
    return (1 - (alpha*p)**(n-1)/(1+(alpha*p)**n)**m)**2/(1+(alpha*p)**n)**(m/2)

# pressure range
p = np.linspace(0,5,100) 

# plot
fig = plt.figure()
ax = fig.add_subplot()

for x in range(nb):
    alpha = alpha_val[x]
    n = n_val[x]
    ax.plot(p,kr(p,alpha,n),label=f'(alpha,n)=({alpha},{n})')

ax.set_xlabel('p')
ax.set_ylabel('kr')
ax.legend()

#plt.show()
plt.savefig(f'{path_fig}relativepermeability.png', dpi = 300)
