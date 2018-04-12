import glob
import matplotlib.pyplot as plt
from os import path
from netCDF4 import Dataset
from scipy.stats import t as tdist
import numpy as np

methodDict = {'ARS232-native':  7,
              'KGS242-native':  8,
              'KGS252-native':  9,
              'KGS254-native': 14,
              'ARS232':        22,
              'DBM453':        23,
              'ARS222':        24,
              'ARS233':        25,
              'ARS343':        26,
              'ARS443':        27,
              'ARK324':        28,
              'ARK436':        29,
              'SSP3333b':      30,
              'SSP3333c':      31}

# Acceptable relative error
errorMax = 0.05

# Load reference solution
fileName = 'output_tsteptype5_tstep30_hydrostatic/nonhydro-X1-dcmip2012_test41.nc'
print "reading reference solution from " + fileName
data = Dataset(fileName)
uRef = data['u'][:]
u0 = np.mean(uRef[0,:,:,:])

# Iterate through methods and find largest step size
maxsizeDict = {}
maxerrorDict = {}
for m,method in enumerate(methodDict.keys()):
  maxsizeDict[method] = -np.inf
  maxerrorDict[method] = np.inf
  print method
  # Load timestep and solution data
  globstr = 'output_tsteptype%s_tstep*/nonhydro-X1-dcmip2012_test41.nc' % methodDict[method]
  for fileName in glob.glob(globstr):
    print "reading " + fileName
    data = Dataset(fileName)
    u = data['u'][:]
    if (np.shape(u)[0] == 31):
      error = np.sqrt( np.sum( (u-uRef)**2, axis=(1,2,3) ) / \
                        (np.shape(u)[1]*np.shape(u)[2]*np.shape(u)[3]) ) / u0
      if (np.amax(error) < errorMax):
        words = fileName.split('_')
        words = words[2].split('/')
        dt = float(words[0].replace('tstep',''))
        if (dt > maxsizeDict[method]):
          maxsizeDict[method] =  dt
          maxerrorDict[method] = np.amax(error)
      else:
        print '... skipping due to incomplete results ...'

print ""
for method in methodDict.keys():
  print method, methodDict[method], maxsizeDict[method], maxerrorDict[method]
