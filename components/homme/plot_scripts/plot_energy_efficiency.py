import glob
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import numpy as np
import os

matplotlib.rcParams.update({'font.size':22})

methodDict = {#'KGU35-native': 5,
              #'ARS232-native':  7,
              #'KGS252-native':  9,
              #'KGS254-native': 14,
              'KGU35': 21,
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

X = 1


# Iterate through methods and plot energy efficiency test
f, ax = pyplot.subplots(figsize=(10,10))
f.set_tight_layout(True)
for m,method in enumerate(methodDict.keys()):
  relErrorDict = {}
  walltimeDict = {}
  print method
  globstr = 'tsteptype%d_tstep*_X%d.out' % (methodDict[method], X)
  for fileName in glob.glob(globstr):
    words = fileName.split('_')
    dt = words[1].replace('tstep','')
    dt = dt.replace('.out','')
    print 'Reading diagnostics in ' + fileName
    fileObj = open(fileName)
    lines = list(fileObj)
    fileObj.close()
    flag = False
    for line in reversed(lines):
      if ('Finished main timestepping loop' in line):
        flag = True
      if (flag and '(E-E0)/E0' in line):
        words = line.split()
        relError = float(words[1])
        relErrorDict[dt] = relError
        break
    if (not flag):
      print '... skipping due to incomplete results ...'
      continue
    fileObj = open(fileName.replace('.out','.err'))
    lines = list(fileObj)
    fileObj.close()
    words = lines[1].split()
    word = words[1] # ignore the word 'real'
    words = word.split(('m'))
    seconds = int(words[0])*60 + float(words[1].replace('s',''))
    walltimeDict[dt] = seconds

  # if no solutions were found, goto next method
  dtList = relErrorDict.keys()
  if (len(dtList) == 0):
    continue

  # need dictionary to map strings to floats
  dtDict = {}
  for dt in dtList:
    dtDict[float(dt)] = dt

  # Sort values and add to plot
  dtPlot = np.sort(dtDict.keys())
  relErrorPlot = np.empty(len(dtPlot))
  walltimePlot = np.empty(len(dtPlot))
  for j,dt in enumerate(dtPlot):
    relErrorPlot[j] = abs(relErrorDict[dtDict[dt]])
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

  ax.loglog(walltimePlot, relErrorPlot, lineStyle, label=method, linewidth=2, markersize=12)

ax.set_ylabel('Relative Energy Error', fontsize='xx-large')
ax.set_xlabel('wall time (s)', fontsize='xx-large')

f.savefig('efficiencyEnergy.png')

ax.legend(loc='best')
f.savefig('efficiencyEnergy_legend.png')
pyplot.show()
