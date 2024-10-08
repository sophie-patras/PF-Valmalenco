"""
For PF postprocessing
Pressure evolution

Created on Wed Aug 7 16:01:20 2024
@author: S.P.
"""

# import os
import numpy as np
import matplotlib.pyplot as plt


from parflow.tools.fs import get_absolute_path
from parflow.tools.io import read_pfb


satur_data = read_pfb(get_absolute_path("../RunProcess/Tmp/PLT.out.satur.00010.pfb"))

# print(f'Dimensions of output file: {press_data.shape}') # plot (NZ,NY,NX)
# print(type(press_data))  # numpy array

satur_shape = satur_data.shape

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x,y = np.meshgrid(np.arange(saturshape[2]), np.arange(saturshape[1]))

for p in range(saturshape[0]):
	saturation = satur_data[p]
	ax.plot_surface(x,y,np.full_like(saturation, p), facecolors=plt.cm.viridis(saturation), rstride=1, cstride=1, antialiased=True, shade=False)

print(saturation)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('s')
ax.set_title('saturation')

m = plt.cm.ScalarMappable(cmap=plt.cm.viridis) #, norm=surf.norm)
vmin = N.min()
vmax = N.max()
m.set_clim(vmin,vmax)
plt.colorbar(m, ax=plt.gca())
# plt.legend('sat')
plt.show()
# plt.savefig('/mnt/c/Users/User/Documents/POLIMI/0_TESI/')
