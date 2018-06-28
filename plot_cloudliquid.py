import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np


cmap = 'Greys'
configList = ['RKZ_P3_lmt4_ne30_ne30_FC5AQUAP_intel_cori-knl', \
              'RKZ_P3_lmt4_ne30_ne30_FC5AQUAP_intel_cori-knl']

f, ax = pyplot.subplots(len(configList), 4, sharex=True, sharey=True, figsize=(15,10))
for j, config in enumerate(configList):
  fileName = './%s/DT0001_01_dycore_mac/DT0001_01_dycore_mac.cam.h0.0001-01-01-00000.nc' % config
  print('Reading in ' + fileName)
  data = Dataset(fileName)
  lat = data['lat'][:]
  lon = data['lon'][:]
  ql = data['CLDLIQ'][:]
  qlTotalIC = np.sum(ql, axis=1)
  qlTotalMax = np.amax(qlTotalIC)
  fileName = './%s/DT0001_01_dycore_mac/DT0001_01_dycore_mac.cam.h0.0001-01-01-03600.nc' % config
  print('Reading in ' + fileName)
  data = Dataset(fileName)
  ql = data['CLDLIQ'][:]
  qlTotal = np.sum(ql, axis=1)
  qlTotalMax = max(qlTotalMax, np.amax(qlTotal))
 
  tmp = ax[j,0].scatter(lon, lat, c=qlTotalIC[0,:], cmap=cmap, vmin=0.0, vmax=qlTotalMax)
  ax[j,0].set_title('Initial Condition')
  ax[j,0].set_ylabel(config)
  ax[j,1].scatter(lon, lat, c=qlTotal[0,:], cmap=cmap, vmin=0.0, vmax=qlTotalMax)
  ax[j,1].set_title('dt=1s')

  for k, dt in enumerate([30, 1800]):
    fileName = './%s/DT%04d_01_dycore_mac/DT%04d_01_dycore_mac.cam.h0.0001-01-01-03600.nc' % (config, dt, dt)
    print('Reading in ' + fileName)
    data = Dataset(fileName)
    ql = data['CLDLIQ'][:]
    ax[j,k+2].scatter(lon, lat, c=qlTotal[0,:], cmap=cmap, vmin=0.0, vmax=qlTotalMax)
    ax[j,k+2].set_title('dt=%ss' % dt)

f.subplots_adjust(right=0.8)
cax = f.add_axes([0.85, 0.15, 0.05, 0.7])
pyplot.colorbar(tmp,cax=cax)
pyplot.show()  
  

