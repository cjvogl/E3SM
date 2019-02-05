import matplotlib
import matplotlib.pyplot as pyplot
from mpl_toolkits.basemap import Basemap
import numpy as np


f1, ax1 = pyplot.subplots(figsize=(10,6))
f2, ax2 = pyplot.subplots(figsize=(10,6))
f3, ax3 = pyplot.subplots(figsize=(10,6))
bmap = Basemap(lon_0=180)

tmp = np.loadtxt('./data/surface_pressure_reference_data.txt')
X = tmp[0:256:]
Y = tmp[256:512,:]
val = tmp[512:768,:]

bmap.ax = ax1
tmp = bmap.pcolor(X, Y, val, cmap='bwr', vmin=-10, vmax=10)
cb = f1.colorbar(tmp, ax=ax1)
bmap.drawcoastlines()

tmp = np.loadtxt('./data/surface_pressure_tolerance_spatial_data.txt')
X = tmp[0:256:]
Y = tmp[256:512,:]
val = tmp[512:768,:]

bmap.ax = ax2
tmp = bmap.pcolor(X, Y, val, cmap='bwr', vmin=-0.05, vmax=0.05)
f2.colorbar(tmp, ax=ax2)
bmap.drawcoastlines()

tmp = np.loadtxt('./data/surface_pressure_tolerance_rms_data.txt')
t = tmp[0,:]
val = tmp[1,:]

ax3.plot(t, val, '-o', color='black', linewidth=3)
ax3.set_xlabel('time (days)', fontsize='xx-large')
ax3.set_ylabel('RMS tolerance (kPa)', fontsize='xx-large')
ax3.set_xlim(0,15)

f1.tight_layout()
f2.tight_layout()
f3.tight_layout()

f1.savefig('surface_pressure_reference.pdf')
f2.savefig('surface_pressure_tolerance_spatial.pdf')
f3.savefig('surface_pressure_tolerance_rms.pdf')
pyplot.show()
