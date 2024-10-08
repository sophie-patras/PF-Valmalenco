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

# plotParams
# text_kwargs = dict(ha='center', va='center', fontsize=12, color='k')

fn = "/home/patras/Valmalenco/Tmp/PLT.out.mask.pfb"
npdata = read_pfb(get_absolute_path(fn))

datashape = npdata.shape
print(f'Dimensions of output file: {datashape}') # plot (NZ,NY,NX)

x,y = np.meshgrid(np.arange(datashape[2]), np.arange(datashape[1]))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for p in range(datashape[0]):
    data = npdata[p]
    print(data)
    ax.plot_surface(x,y,np.full_like(data, p), facecolors=plt.cm.viridis(data), rstride=1, cstride=1, antialiased=True, shade=False)

ax.set_xlabel('x')
ax.set_ylabel('y')
#ax.set_zlabel('p')
ax.set_title('mask')

m = plt.cm.ScalarMappable(cmap=plt.cm.viridis)  #, norm=surf.norm)
vmin = npdata.min()
vmax = npdata.max()
m.set_clim(vmin,vmax)
plt.colorbar(m, aspect = 10, pad=0.1, fraction=0.05, ax=plt.gca())

# ax.text2D(0.1,0.9,f"time={DumpInt*DumpGap}h",transform = ax.transAxes)
plt.show()
#plt.pause(0.5)

# plt.savefig('/mnt/c/Users/User/Documents/POLIMI/0_TESI/')
