from dictionaries import create_walltime_dict, \
                         methodDict, lineStyleDict, colorDict
import matplotlib.pyplot as pyplot
import matplotlib
from netCDF4 import Dataset
import numpy as np
import os

matplotlib.rcParams.update({'font.size':22})

dtBuffer = 10
timeBuffer = 100

# Iterate through methods and plot energy efficiency test
fSize, axSize = pyplot.subplots(figsize=(10,5))
fStageSize, axStageSize = pyplot.subplots(figsize=(10,5))
fTime, axTime = pyplot.subplots(figsize=(10,5))
fShortestTime, axShortestTime = pyplot.subplots(figsize=(10,5))
fSize.set_tight_layout(True)
fStageSize.set_tight_layout(True)
fTime.set_tight_layout(True)
fShortestTime.set_tight_layout(True)
count = 0
largestDtDict = {}
largestStageDtDict = {}
timeDict = {}
shortestTimeDict = {}
for method in methodDict.keys():
  print method
  globstr = 'tsteptype%d_tstep*.out' % methodDict[method]
  walltimeDict = create_walltime_dict(globstr, 0, np.inf, '', 'restart')

  # if no solutions were found, goto next method
  dtList = walltimeDict.keys()
  if (len(dtList) == 0):
    continue

  # need dictionary to map strings to floats
  dtDict = {}
  for dt in dtList:
    dtDict[float(dt)] = dt

  dt = np.amax(np.sort(dtDict.keys()))
  largestDtDict[method] = dt
  if ('o' in lineStyleDict[method]):
    iStageNum = 1
  elif ('x' in lineStyleDict[method]):
    iStageNum = 2
  elif ('^' in lineStyleDict[method]):
    iStageNum = 3
  elif ('s' in lineStyleDict[method]):
    iStageNum = 4
  elif ('p' in lineStyleDict[method]):
    iStageNum = 5
  else:
    exit('error')
  largestStageDtDict[method] = dt/iStageNum
  if ('KGU35' in method):
    timeDict[method] = 0.0
  else:
    timeDict[method] = walltimeDict['100']
  shortestTimeDict[method] = walltimeDict[dtDict[dt]]

# bubble sort largestDtDict and then plot
methodList = largestDtDict.keys()
dtList = largestDtDict.values()
for m in range(len(dtList)):
  for n in range(0,len(dtList)-m-1):
    if (dtList[n] > dtList[n+1]):
      swap = dtList[n]
      dtList[n] = dtList[n+1]
      dtList[n+1] = swap
      swap = methodList[n]
      methodList[n] = methodList[n+1]
      methodList[n+1] = swap
for m in range(len(dtList)):
  dt = dtList[m]
  method = methodList[m]
  axSize.plot((0.0, dt), (m, m), lineStyleDict[method], label=method, \
              color=colorDict[method], linewidth=3, markersize=12)
  axSize.text(dt+dtBuffer, m, '%d' % dt, horizontalalignment='left',
                         verticalalignment='center', color='k',
                         clip_on=True)

# bubble sort largestStageDtDict and then plot
methodList = largestStageDtDict.keys()
dtList = largestStageDtDict.values()
for m in range(len(dtList)):
  for n in range(0,len(dtList)-m-1):
    if (dtList[n] > dtList[n+1]):
      swap = dtList[n]
      dtList[n] = dtList[n+1]
      dtList[n+1] = swap
      swap = methodList[n]
      methodList[n] = methodList[n+1]
      methodList[n+1] = swap
for m in range(len(dtList)):
  dt = dtList[m]
  method = methodList[m]
  axStageSize.plot((0.0, dt), (m, m), lineStyleDict[method], label=method, \
              color=colorDict[method], linewidth=3, markersize=12)
  axStageSize.text(dt+dtBuffer, m, '%d' % dt, horizontalalignment='left',
                         verticalalignment='center', color='k',
                         clip_on=True)

# bubble sort timeDict and then plot
methodList = timeDict.keys()
timeList = timeDict.values()
for m in range(len(timeList)):
  for n in range(0,len(timeList)-m-1):
    if (timeList[n] < timeList[n+1]):
      swap = timeList[n]
      timeList[n] = timeList[n+1]
      timeList[n+1] = swap
      swap = methodList[n]
      methodList[n] = methodList[n+1]
      methodList[n+1] = swap
for m in range(len(dtList)):
  time = timeList[m]
  method = methodList[m]
  axTime.plot((0.0, time), (m, m), lineStyleDict[method], label=method, \
            color=colorDict[method], linewidth=3, markersize=12)
  axTime.text(time+timeBuffer, m, '%d' % time, horizontalalignment='left',
                         verticalalignment='center', color='k',
                         clip_on=True)


# bubble sort shortestTimeDict and then plot
methodList = shortestTimeDict.keys()
timeList = shortestTimeDict.values()
for m in range(len(timeList)):
  for n in range(0,len(timeList)-m-1):
    if (timeList[n] < timeList[n+1]):
      swap = timeList[n]
      timeList[n] = timeList[n+1]
      timeList[n+1] = swap
      swap = methodList[n]
      methodList[n] = methodList[n+1]
      methodList[n+1] = swap
for m in range(len(dtList)):
  time = timeList[m]
  method = methodList[m]
  axShortestTime.plot((0.0, time), (m, m), lineStyleDict[method], label=method, \
            color=colorDict[method], linewidth=3, markersize=12)
  axShortestTime.text(time+timeBuffer, m, '%d' % time, horizontalalignment='left',
                         verticalalignment='center', color='k',
                         clip_on=True)

axSize.set_xlabel('timestep (s)', fontsize='xx-large')
axStageSize.set_xlabel('timestep (s)', fontsize='xx-large')
axTime.set_xlabel('walltime for dt=100s (s)', fontsize='xx-large')
axShortestTime.set_xlabel('shortest walltime (s)', fontsize='xx-large')

axSize.legend(loc='best')
axStageSize.legend(loc='best')
axTime.legend(loc='best')
axShortestTime.legend(loc='best')
pyplot.show()
