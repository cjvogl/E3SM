from dictionaries import methodDict
from netCDF4 import Dataset
import numpy as np
import sys
import os

if (len(sys.argv) < 2):
  print "Usage: python calc_largest_coupling_timestep_X10_baroclinicinstability.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

methodList = ['ARS232', 'ARS222', 'ARS233', 'ARS343',
              'ARS443', 'ARK324', 'ARK436', 'DBM453', 'SSP3333b', 'SSP3333c',
              'IMKG232a', 'IMKG232b', 'IMKG242a', 'IMKG242b', 'IMKG243a', 
              'IMKG252a', 'IMKG252b', 'IMKG253a', 'IMKG253b', 'IMKG254a', 
              'IMKG254b', 'IMKG254c', 'IMKG343a', 'GSA222', 'SSP2232',
              'ARK548', 'ARK437']

dtList = [300, 270, 240, 216, 200, 192, 180, 160, 150, 135, 120, 100, 50, 20 ,10]

# Iterate through methods and find largest timestep
largestDtDict = {}
for method in methodList:
  print method
  currentList = dtList + [] # make a copy of dtList

  success = False
  dt = np.amax(currentList)
  exceed = np.inf
  while (not success and dt != 0.0):
    if (dt % 10 != 0):
      dtString = '%0.1f' % (dt/10.0)
    else:
      dtString = '%d' % (dt/10.0)
    fileName = outpath + '/output_tsteptype%d_tstep%s_dcmip4_X10/dcmip2012_test41.nc' \
               % (methodDict[method], dtString)
    if (not os.path.exists(fileName)):
      data = {'time': [0.0,]}
    else:
      data = Dataset(fileName)
    print('Checking if ' + fileName + ' completed...')
    if (data['time'][-1] < 1.45):
      print('... no.')
      for j, dtau in enumerate(currentList):
        if (dtau == dt):
          currentList[j] = 0.0
          dt = np.amax(currentList)
          break
    else:
      print('... yes.')
      success = True

  largestDtDict[method] = dt/10.0

filename = './data/largest_coupling_timestep_X10_baroclinicinstability.txt'
print('Writing out ' + filename)
fileobj = open(filename, 'w')
for method in methodList:
  fileobj.write('%11s %0.1f\n' % (method, largestDtDict[method]))
fileobj.close()
