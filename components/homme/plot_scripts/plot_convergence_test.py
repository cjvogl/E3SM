import glob
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import numpy as np
from dictionaries import methodDict, lineStyleDict, colorDict

matplotlib.rcParams.update({'font.size':22})

# If you don't want all available results to be plotted, specify which methods# to plot here

# methodDict = {'KGU35-native': 5}

def plot_convergence_test(fileRef, indRef, varRef, \
                          suffix='', suffix_omit=None, minDt=0.0, maxDt=np.inf):
  # Parse info from reference solution format: ./tsteptype*_tstep(dtRef)*/(testName)
  words = fileRef.split('/')
  testName = words[-1]
  words = fileRef.split('_')
  i = 0
  for word in words:
    i = i+1
    if ('output' in word):
      break
  dtRef = float(words[i+1].replace('tstep',''))

  # Read in reference solution
  print 'Reading reference solution from ' + fileRef
  data = Dataset(fileRef)
  qRef = data[varRef][:]
  qRef = qRef[indRef,:,:,:]
  tRef = data['time'][:]
  tRef = tRef[indRef]
  shape = np.shape(qRef)
  numElements = shape[0]*shape[1]*shape[2]

  # Create figures for RMS and Linf plots
  f1, ax1 = pyplot.subplots(figsize=(10,10))
  f2, ax2 = pyplot.subplots(figsize=(10,10))

  # Iterate through methods
  for m,method in enumerate(methodDict.keys()):
    solutionDict = {}
    print method
    # Parse all stdout files for current method and given suffix
    globstr = 'tsteptype%d_tstep*%s.out' % (methodDict[method], suffix)
    for fileName in glob.glob(globstr):
      if (suffix_omit is not None and suffix_omit in fileName):
        continue
      words = fileName.split('_')
      dt = words[1].replace('tstep','')
      dt = dt.replace('.out','')
      # If current timestep is appropriate, obtain solution if it exists
      if (float(dt) > dtRef+1e-12 and float(dt) < maxDt and float(dt) > minDt):
        directory = './output_'+fileName.replace('.out','')
        print 'Reading solution in ' + directory
        data = Dataset(directory+'/'+testName)
        q = data[varRef][:]
        t = data['time'][:]
        if (len(t) > indRef and abs(t[indRef] - tRef) < 1e-10):
          solutionDict[dt] = q[indRef,:,:,:]
        else:
          print '... skipping due to incomplete results ...'

    # If no solutions were found for current method, goto next method
    dtList = solutionDict.keys()
    if (len(dtList) == 0):
      continue

    # Compute errors from reference solution
    dtDict = {}
    RMSerror = {}
    LIerror = {}
    for dt in dtList:
      q = solutionDict[dt]
      dtDict[float(dt)] = dt
      RMSerror[dt] = np.sqrt(np.sum((q-qRef)**2)/numElements)
      LIerror[dt] = np.amax(abs(q-qRef))
    # Sort values by timestep for plotting
    dtPlot = np.sort(dtDict.keys())
    RMSPlot = np.empty(len(dtPlot))
    LIPlot = np.empty(len(dtPlot))
    for j,dt in enumerate(dtPlot):
      RMSPlot[j] = RMSerror[dtDict[dt]]
      LIPlot[j] = LIerror[dtDict[dt]]
    # Approximate order of convergence
    if (len(dtList) > 1):
      orderL2 = np.log(RMSPlot[1:]/RMSPlot[0:-1])/np.log(dtPlot[1:]/dtPlot[0:-1])
      orderLI = np.log(LIPlot[1:]/LIPlot[0:-1])/np.log(dtPlot[1:]/dtPlot[0:-1])
    else:
      orderL2 = (0.0,)
      orderLI = (0.0,)
    # Plot RMS and Linf errors
    ax1.loglog(dtPlot, RMSPlot, lineStyleDict[method], \
               label='%s (final=%3.2f, best=%3.2f)' % (method,orderL2[0],np.amax(orderL2)), \
               color=colorDict[method], linewidth=3, markersize=12)
    ax2.loglog(dtPlot, LIPlot, lineStyleDict[method], \
               label='%s (final=%3.2f, best=%3.2f)' % (method,orderLI[0],np.amax(orderLI)), \
               color=colorDict[method], linewidth=3, markersize=12)
  # end of method loop

  # Finalize figures
  ax1.set_ylabel('RMS error', fontsize='xx-large')
  ax1.set_xlabel('dt (s)', fontsize='xx-large')
  ax1.axis('equal')
  ax2.set_ylabel('LI Error (%s)' %varRef, fontsize='xx-large')
  ax2.set_xlabel('dt (s)', fontsize='xx-large')
  ax2.axis('equal')
  f1.tight_layout()
  f2.tight_layout()

  # Save figures without legends
  f1.savefig('errorRMS%s_%s.png' % (varRef, suffix))
  f2.savefig('errorLI%s_%s.png' % (varRef, suffix))

  # Save figures with legends
  ax1.legend(loc='best')
  ax2.legend(loc='best')
  f1.savefig('errorRMS%s_%s_legend.png' % (varRef, suffix))
  f2.savefig('errorLI%s_%s_legend.png' % (varRef, suffix))
  pyplot.show()

if (__name__ == '__main__'):
  import os

  cwd = os.getcwd()
  if ('gravitywave_test' in os.getcwd()):
    testName = 'dcmip2012_test31'
    fileRef = cwd + '/output_tsteptype5_tstep0.000390625_nmax768000_nu0.0/' \
                  + 'dcmip2012_test31.nc'
    indRef = 10
    varRef = 'T'
    suffix = 'nu0.0'
    plot_convergence_test(fileRef, indRef, varRef, suffix=suffix)
    fileRef = cwd + '/output_tsteptype5_tstep0.000390625_nmax768000/' \
                  + 'dcmip2012_test31.nc'
    plot_convergence_test(fileRef, indRef, varRef, suffix_omit=suffix)
    pyplot.show()
