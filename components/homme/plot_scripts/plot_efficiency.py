import glob
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import numpy as np
import os

noHV = False
rtol = '1e-3'
splitting = '1'

matplotlib.rcParams.update({'font.size':22})

methodDict = {'KGU35': 21,
              'ARS232': 22,
              'DBM453': 23,
              'ARS222': 24,
              'ARS233': 25,
              'ARS343': 26,
              'ARS443': 27,
              'ARK324': 28,
              'ARK436': 29,
              'SSP3333b': 30,
              'SSP3333c': 31,
              'KGS232': 32,
              'KGS242': 33,
              'KGS243': 34,
              'KGS252': 35,
              'KGS254': 36}

orderDict = {'KGU35-native': 3,
              'ARS232-native': 2,
              'KGS252-native': 2,
              'KGS254-native': 2,
              'KGU35': 3,
              'ARS232': 2,
              'DBM453': 3,
              'ARS222': 2,
              'ARS233': 3,
              'ARS343': 3,
              'ARS443': 3,
              'ARK324': 3,
              'ARK436': 4,
              'SSP3333b': 3,
              'SSP3333c': 3,
              'KGS232': 2,
              'KGS242': 2,
              'KGS243': 2,
              'KGS252': 2,
              'KGS254': 2}

istageDict = {'KGU35-native': 0,
              'ARS232-native':  2,
              'KGS252-native':  2,
              'KGS254-native': 4,
              'KGU35': 0,
              'ARS232': 2,
              'DBM453': 4,
              'ARS222': 2,
              'ARS233': 2,
              'ARS343': 3,
              'ARS443': 4,
              'ARK324': 3,
              'ARK436': 5,
              'SSP3333b': 2,
              'SSP3333c': 2,
              'KGS232': 2,
              'KGS242': 2,
              'KGS243': 3,
              'KGS252': 2,
              'KGS254': 4}

# Load reference solution
if ('gravitywave_test' in os.getcwd()):
  testName = 'dcmip2012_test31'
  dtRef = 0.000390625
  nmax = 768000
  indRef = 10
  varRef = 'T'
  if (noHV):
    directory = './output_tsteptype5_tstep%10.9f_nmax%d_nu0.0/%s.nc' \
                % (dtRef, nmax, testName)
  else:
    directory = './output_tsteptype5_tstep%10.9f_nmax%d/%s.nc' \
                % (dtRef, nmax, testName)
elif ('baroclinicinstability_test' in os.getcwd()):
  testName = 'dcmip2012_test41'
  dtRef = 10
  indRef = 15
  varRef = 'u'
  if (noHV):
    directory = './output_tsteptype5_tstep%2.0f_hydrostatic_nu0.0/%s.nc' \
                % (dtRef, testName)
  else:
    directory = './output_tsteptype5_tstep%2.0f_hydrostatic/%s.nc' \
                % (dtRef, testName)

print 'Reading reference solution from ' + directory
data = Dataset(directory)
qRef = data[varRef][:]
qRef = qRef[indRef,:,:,:]
tRef = data['time'][:]
tRef = tRef[indRef]
shape = np.shape(qRef)
numElements = shape[0]*shape[1]*shape[2]

# Iterate through methods and plot convergence test
f1, ax1 = pyplot.subplots(figsize=(10,10))
f2, ax2 = pyplot.subplots(figsize=(10,10))
for m,method in enumerate(methodDict.keys()):
  dtList = []
  solutionDict = {}
  walltimeDict = {}
  print method
  globstr = 'tsteptype%d_tstep*.out' % methodDict[method]
  for fileName in glob.glob(globstr):
    if ((noHV and "nu0.0" in fileName) or (not noHV and "nu0.0" not in fileName)):
      words = fileName.split('_')
      dt = words[1].replace('tstep','')
      dt = dt.replace('.out','')
      if (float(dt) > dtRef+1e-12):
        directory = './output_'+fileName.replace('.out','')
        print 'Reading solution in ' + directory
        data = Dataset(directory+'/'+testName+'.nc')
        q = data[varRef][:]
        t = data['time'][:]
        if (len(t) > indRef and abs(t[indRef] - tRef) < 1e-10):
          dtList.append(dt)
          solutionDict[dt] = q[indRef,:,:,:]
        else:
          print '... skipping due to incomplete results ...'
          continue
        f = open(fileName.replace('.out','.err'))
        lines = list(f)
        f.close()
        words = lines[1].split()
        word = words[1] # ignore the word 'real'
        words = word.split(('m'))
        seconds = int(words[0])*60 + float(words[1].replace('s',''))
        walltimeDict[dt] = seconds

  # if no solutions were found, goto next method
  if (len(dtList) == 0):
    continue

  # Compute errors from reference solution
  dtDict = {}
  L2error = {}
  LIerror = {}
  for dt in dtList:
    q = solutionDict[dt]
    dtDict[float(dt)] = dt
    L2error[dt] = np.sqrt(np.sum((q-qRef)**2)/numElements)
    LIerror[dt] = np.amax(abs(q-qRef))

  # Sort values and add to plot
  dtPlot = np.sort(dtDict.keys())
  L2Plot = np.empty(len(dtPlot))
  LIPlot = np.empty(len(dtPlot))
  walltimePlot = np.empty(len(dtPlot))
  for j,dt in enumerate(dtPlot):
    L2Plot[j] = L2error[dtDict[dt]]
    LIPlot[j] = LIerror[dtDict[dt]]
    walltimePlot[j] = walltimeDict[dtDict[dt]]
  if (orderDict[method] == 2):
    lineStyle = '-'
  elif (orderDict[method] == 3):
    lineStyle = '--'
  elif (orderDict[method] == 4):
    lineStyle = '-.'
  else:
    print '... need to add linestyle for order %d ...' % orderDict[method]
    exit()

  if (istageDict[method] == 0):
    lineStyle = lineStyle + 'o'
  elif (istageDict[method] == 2):
    lineStyle =  lineStyle +'x'
  elif (istageDict[method] == 3):
    lineStyle = lineStyle + '^'
  elif (istageDict[method] == 4):
    lineStyle = lineStyle + 's'
  elif (istageDict[method] == 5):
    lineStyle = lineStyle +'p'
  elif (istageDict[method] == 6):
    lineStyle = lineStyle + 'h'
  else:
    print '... need to add linestyle for istage %s ...' % istageDict[method]
    exit()

  ax1.loglog(walltimePlot, L2Plot, lineStyle, label=method, linewidth=2, markersize=12)
  ax2.loglog(walltimePlot, LIPlot, lineStyle, label=method, linewidth=2, markersize=12)

ax1.set_ylabel('L2 Error', fontsize='xx-large')
ax1.set_xlabel('wall time (s)', fontsize='xx-large')
ax2.set_ylabel('LI Error', fontsize='xx-large')
ax2.set_xlabel('wall time (s)', fontsize='xx-large')

f1.tight_layout()
f1.savefig('efficiencyL2.png')
f2.tight_layout()
f2.savefig('efficiencyLI.png')

ax1.legend(loc='best')
ax2.legend(loc='best')
f1.savefig('efficiencyL2_legend.png')
f2.savefig('efficiencyLI_legend.png')
pyplot.show()
