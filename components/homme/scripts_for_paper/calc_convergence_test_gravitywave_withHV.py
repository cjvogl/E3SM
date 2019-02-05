import glob
from dictionaries import methodDict, lineStyleDict, colorDict
from netCDF4 import Dataset
import numpy as np
import sys

if (len(sys.argv) < 2):
  print "Usage: python calc_convergence_test_gravitywave_withHV.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

methodList = ['KGU35', 'ARS232', 'DBM453', 'ARS222', 'ARS233', 'ARS343', \
              'ARS443', 'ARK324', 'ARK436', 'SSP3333b', 'SSP3333c']

# read in reference solution
fileRef = outpath + '/output_tsteptype5_tstep0.000390625_nmax768000/dcmip2012_test31.nc'
print 'Reading reference solution from ' + fileRef
data = Dataset(fileRef)
qRef = data['T'][:]
qRef = qRef[10,:,:,:]
tRef = data['time'][:]
tRef = tRef[10]
# iterate through methods and plot convergence test
for method in methodList:
  solutionDict = {}
  print method
  # iterate through all .out files
  globstr = outpath + '/tsteptype%d_tstep*.out' % methodDict[method]
  for filePath in glob.glob(globstr):
    # skip runs without viscosity
    if ('nu0.0' in filePath):
      continue
    # obtain the timestep
    words = filePath.split('/')
    fileName = words[-1]
    words = fileName.split('_')
    dt = words[1].replace('tstep','')
    # check that timestep is small enough to be in the asymptotic regime
    if (float(dt) < 1.1):
      directory = outpath + '/output_'+fileName.replace('.out','')
      print 'Reading solution in ' + directory
      data = Dataset(directory+'/dcmip2012_test31.nc')
      q = data['T'][:]
      t = data['time'][:]
      solutionDict[dt] = q[10,:,:,:]

  # compute maximum relative error from reference solution
  dtList = solutionDict.keys()
  dtDict = {}
  error = {}
  for dt in dtList:
    q = solutionDict[dt]
    dtDict[float(dt)] = dt
    error[dt] = np.amax(abs(q-qRef)/qRef)
  # sort values by timestep for plotting
  dtPlot = np.sort(dtDict.keys())
  errorPlot = np.empty(len(dtPlot))
  for j,dt in enumerate(dtPlot):
    errorPlot[j] = error[dtDict[dt]]
  output = np.vstack((dtPlot, errorPlot))
  filename = './data/' + method +'_convergence_data_withHV.txt'
  print('Writing out ' + filename)
  np.savetxt(filename, output)
