import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np
import os

wrkdir = os.getenv('S') + '/homme/jw_baroclinic'

data = Dataset(wrkdir + '/movies/jw_baroclinic1.nc',mode='r')
dataARK = Dataset(wrkdir + '/U35_ARK/jw_baroclinic1.nc',mode='r')
var = data.variables['u'][:]
varARK = dataARK.variables['u'][:]

times = (1,5,10)

f,axarray = pyplot.subplots(len(times),3,figsize=(16,4*len(times)))
for i,time in enumerate(times):
    level = np.argmax(abs(var[time,:,:,:] - varARK[time,:,:,:]))
    (level, tmp1, tmp2) = np.unravel_index(level,np.shape(var[time,:,:,:]))
    pvar = var[time,level,:,:]
    pvarARK = varARK[time,level,:,:]

    vmin = np.amin(pvar)
    vmax = np.amax(pvar)
    axarray[i,0].set_title('Original after %d step(s)' % time)
    axarray[i,0].set_ylabel('Level %d' % level)
    cax = axarray[i,0].pcolor(pvar,vmin=vmin,vmax=vmax)
    f.colorbar(cax,ax=axarray[i,0])
    axarray[i,1].set_title('ARKode after %d step(s)' % time)
    cax = axarray[i,1].pcolor(pvarARK,vmin=vmin,vmax=vmax)
    f.colorbar(cax,ax=axarray[i,1])
    axarray[i,2].set_title('Difference after %d step(s)' % time)
    cax = axarray[i,2].pcolor(pvar - pvarARK)
    f.colorbar(cax,ax=axarray[i,2])




pyplot.show()
