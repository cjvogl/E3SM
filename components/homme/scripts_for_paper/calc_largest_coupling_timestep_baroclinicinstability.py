from dictionaries import methodDict
from netCDF4 import Dataset
import numpy as np
import sys

if (len(sys.argv) < 2):
  print "Usage: python calc_largest_coupling_timestep_baroclinicinstability.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

methodList = ['ARS232', 'ARS222', 'ARS233', 'ARS343',
              'ARS443', 'ARK324', 'ARK436', 'DBM453', 'SSP3333b',
              'IMEX-KG232a', 'IMEX-KG232b', 'IMEX-KG242a', 'IMEX-KG242b',
              'IMEX-KG243b', 'IMEX-KG252a', 'IMEX-KG252b',
              'IMEX-KG253a', 'IMEX-KG254a', 'IMEX-KG254b',
              'IMEX-KG343a']

dtList = [300, 270, 240, 216, 200, 192, 180, 160, 150, 135, 120, 100, 50, 20 ,10]

# Iterate through methods and find acceptable timestep
largestDtDict = {}
dtDict = {}
for method in methodList:
  print method
  currentList = dtList + [] # make a copy of dtList

  success = False
  dt = np.amax(currentList)
  while (not success):
    fileName = outpath + '/output_tsteptype%d_tstep%d/dcmip2012_test41.nc' \
               % (methodDict[method], dt)
    data = Dataset(fileName)
    print('Checking if ' + fileName + ' completed...')
    if (data['time'][-1] < 14.5):
      print('... no.')
      for j, dtau in enumerate(currentList):
        if (dtau == dt):
          currentList[j] = 0.0
          dt = np.amax(currentList)
          break
    else:
      print('... yes.')
      success = True

  largestDtDict[method] = dt

filename = './data/largest_coupling_timestep_baroclinicinstability.txt'
print('Writing out ' + filename)
fileobj = open(filename, 'w')
for method in methodList:
  fileobj.write('%s\t%d\n' % (method, largestDtDict[method]))
fileobj.close()
