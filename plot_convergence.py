import glob
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np


#dtList = [1800,450,120,30,8,2,1]
dtList = [1800,450,120,30,8,1]
configDict = {}
resolution = 'ne30_ne30'
compset = 'FC5AQUAP'
compiler = 'intel'
machine = 'edison'



f1, ax1 = pyplot.subplots(figsize=(10,10))
f2, ax2 = pyplot.subplots(figsize=(10,10))
for dt in dtList:
  globstr = 'RKZ_P*_lmt4_ne30_ne30*/DT%04d_01_dycore_mac/*.cam.h0.0001-01-01-03600.nc' % dt
  for fileName in glob.glob(globstr):
    solutionDict = {}
    words = fileName.split('_')
    config = words[0]
    for j in range(1,len(words)):
      if (words[j] == compset):
        break
      config = config + '_' + words[j]
    if (config not in configDict.keys()):
      solutionDict = {}
      configDict[config] = solutionDict
    solutionDict = configDict[config]
    print('Reading solution in ' + fileName)
    data = Dataset(fileName)
    T = data['T'][:]
    solutionDict[dt] = T[0,:,:]
 
dtRef = min(dtList)
for config in configDict.keys(): 
  dtDict = {}
  L2error = {}
  LIerror = {}
  solutionDict = configDict[config]
  qRef = solutionDict[dtRef]
  numElements = np.shape(qRef)[0]*np.shape(qRef)[1]
  for dt in dtList:
    if (dt is dtRef):
      continue
    q = solutionDict[dt]
    dtDict[float(dt)] = dt
    L2error[dt] = np.sqrt(np.sum((q-qRef)**2)/numElements)
    LIerror[dt] = np.amax(abs(q-qRef))
  
  dtPlot = np.sort(dtDict.keys())
  L2Plot = np.empty(len(dtPlot))
  LIPlot = np.empty(len(dtPlot))
  for j,dt in enumerate(dtPlot):
    L2Plot[j] = L2error[dtDict[dt]]
    LIPlot[j] = LIerror[dtDict[dt]]
  
  if (len(dtList) > 1):
    orderL2 = np.log(L2Plot[1:]/L2Plot[0:-1])/np.log(dtPlot[1:]/dtPlot[0:-1])
    orderLI = np.log(LIPlot[1:]/LIPlot[0:-1])/np.log(dtPlot[1:]/dtPlot[0:-1])
  else:
    orderL2 = (0.0,)
    orderLI = (0.0,)
  
  ax1.loglog(dtPlot, L2Plot, '-o', label='%s final=%3.2f, best=%3.2f' % (config, orderL2[0], np.amax(orderL2)), linewidth=3, markersize=12)
  ax2.loglog(dtPlot, LIPlot, '-o', label='%s final=%3.2f, best=%3.2f' % (config, orderLI[0], np.amax(orderLI)), linewidth=3, markersize=12)

ax1.set_ylabel('L2 Error', fontsize='xx-large')
ax1.set_xlabel('dt', fontsize='xx-large')
ax1.axis('equal')

ax2.set_ylabel('LI Error', fontsize='xx-large')
ax2.set_xlabel('dt', fontsize='xx-large')
ax2.axis('equal')

ax1.legend(loc='best', fontsize='xx-large')
ax2.legend(loc='best', fontsize='xx-large')

pyplot.show()






