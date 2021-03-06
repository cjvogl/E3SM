from dictionaries import methodDict
from netCDF4 import Dataset
import numpy as np
import sys

if (len(sys.argv) < 2):
  print "Usage: python calc_largest_coupling_timestep_baroclinicinstability.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

methodList = ['ARS232', 'ARS222', 'ARS233', 'ARS343',
              'ARS443', 'ARK324', 'ARK436', 'DBM453', 'SSP3333b', 'SSP3333c',
              'IMKG232a', 'IMKG232b', 'IMKG242a', 'IMKG242b', 'IMKG243a', 
              'IMKG252a', 'IMKG252b', 'IMKG253a', 'IMKG253b', 'IMKG254a', 
              'IMKG254b', 'IMKG254c', 'IMKG343a', 'GSA222', 'SSP2232',
              'ARK548', 'ARK437']
methodList = ['ARK548',]

dtList = [300, 270, 240, 216, 200, 192, 180, 160, 150, 135, 120, 100, 50, 20 ,10]

# read in reference solution for surface pressure
fileRef = outpath + '/output_tsteptype5_tstep10_hydrostatic/dcmip2012_test41.nc'
print 'Reading reference solution from ' + fileRef
data = Dataset(fileRef)
qRef = data['ps'][:]

# read in acceptable solution for surface pressure
fileRef = outpath + '/output_tsteptype5_tstep300_hydrostatic/dcmip2012_test41.nc'
print 'Reading acceptable solution from ' + fileRef
data = Dataset(fileRef)
q = data['ps'][:]
numElements = np.shape(q)[1]*np.shape(q)[2]
errorTol = np.sqrt( np.sum((q-qRef)**2, axis=(1,2))/numElements )


# Iterate through methods and find largest timestep
largestDtDict = {}
exceedDict = {}
for method in methodList:
  print method
  currentList = dtList + [] # make a copy of dtList

  if (method == 'ARK548'):
    currentList.append(320)

  success = False
  dt = np.amax(currentList)
  exceed = np.inf
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
      q = data['ps'][:]
      error = np.sqrt( np.sum((q-qRef)**2, axis=(1,2))/numElements )
      ind = np.argmax(error-errorTol)
      exceed = (error[ind]-errorTol[ind])/(error[ind]+1e-16)
      success = True

  largestDtDict[method] = dt
  exceedDict[method] = exceed

filename = './data/largest_coupling_timestep_baroclinicinstability.txt'
print('Writing out ' + filename)
fileobj = open(filename, 'w')
for method in methodList:
  fileobj.write('%11s %3d %2.0f\n' % (method, largestDtDict[method], exceedDict[method]*100))
fileobj.close()
