"""
Relative Permeability, as a function of pressure
Van Genuchten specification, Ref [VG80]
cf PF-UM

Created on Aug 28 2024, 3:41pm
author : SP.
"""

import matplotlib.pyplot as plt
import numpy as np

alpha_val = 1.0
n_val = 1.4

# VanGenuchten
def kr(p,alpha,n):
    m = 1-1/n
    return (1 - (alpha*p)**(n-1)/(1+(alpha*p)**n)**m)**2/(1+(alpha*p)**n)**(m/2)

# pressure range
p = np.linspace(0,5,100) 

# plot
fig = plt.figure()
ax = fig.add_subplot()

ax.plot(p,kr(p,alpha_val,n_val))
ax.set_xlabel('p')
ax.set_ylabel('kr')

plt.show()
