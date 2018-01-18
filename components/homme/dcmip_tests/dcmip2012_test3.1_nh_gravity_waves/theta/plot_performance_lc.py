import glob
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np

noHV = False
rtol = '1e-3'
splitting = '1'

methodDict = {#'U35-ref': 5,
                  #'U35': 11,
                  'ARS232-ref': 7,
                  'ARS232': 12}
                  #'DBM453': 13}
                  #'ARS222': 14,
                  #'ARS233': 15,
                  #'ARS343': 16,
                  #'ARS443': 17,
                  #'ARK324': 18}
                  #'ARK436': 19}

# Iterate through methods and plot nonlinear iteration count
f, ax1 = pyplot.subplots(figsize=(10,10))
f, ax2 = pyplot.subplots(figsize=(10,10))
maxPlotMax = 0
avgPlotMax = 0.0
for m,method in enumerate(methodDict.keys()):
  dtList = []
  maxCountDict = {}
  avgCountDict = {}
  print method
  # Load timestep and solution data
  if (methodDict[method] < 12):
    globstr = 'tsteptype%d_tstep*.out' % methodDict[method]
  else:
    globstr = 'tsteptype%d_tstep*_rtol%s*_splitting%s*.out' % \
                                         (methodDict[method], rtol, splitting)
  for fileName in glob.glob(globstr):
    if ((noHV and "nu0.0" in fileName) or (not noHV and "nu0.0" not in fileName)):
      words = fileName.split('_')
      dt = words[1].replace('tstep','')
      # obtain nonlinear iteration count
      f = open(fileName)
      lines = list(f)
      f.close
      count = 0
      for line in reversed(lines):
        if (count == 2):
          break
        elif (line.startswith('  Max num nonlin iters')):
          words = line.split()
          itercount = int(words[5])
          if (float(dt) not in dtList):
            dtList.append(float(dt))
          maxCountDict[dt] = itercount
          maxPlotMax = max(maxPlotMax,itercount)
          count += 1
        elif (line.startswith('  Avg num nonlin iters')):
          words = line.split()
          itercount = float(words[5])
          if (float(dt) not in dtList):
            dtList.append(float(dt))
          avgCountDict[dt] = itercount
          avgPlotMax = max(avgPlotMax,itercount)
          count += 1

  # Sort values and add to plot
  dtPlot = np.sort(dtList)
  maxCountPlot = np.empty(len(dtPlot))
  avgCountPlot = np.empty(len(dtPlot))
  for j,dt in enumerate(dtPlot):
    maxCountPlot[j] = maxCountDict[str(dt)]
    avgCountPlot[j] = avgCountDict[str(dt)]
  if (m < 7):
    lineStyle = '-o'
  else:
    lineStyle = '-d'
  ax1.semilogx(dtPlot,maxCountPlot,lineStyle,label='%s' % method)
  ax2.semilogx(dtPlot,avgCountPlot,lineStyle,label='%s' % method)


#orderList = ('first', 'second', 'third')
#for j, order in enumerate(orderList):
#    coeff = L2Plot[0]/dtPlot[0]**(j+1)
#    ax1.loglog(dtPlot,coeff*dtPlot**(j+1),'--',label=order)
#    coeff = LIPlot[0]/dtPlot[0]**(j+1)
#    ax2.loglog(dtPlot,coeff*dtPlot**(j+1),'--',label=order)

ax1.set_title('Maximum Nonlinear Iteration Count per Timestep',fontsize=16)
ax1.set_xlabel('dt',fontsize=16)
ax1.set_ylim(0,maxPlotMax+1)
ax1.legend(loc='best',fontsize=16)
ax1.tick_params(axis='both', which='major', labelsize=14)

ax2.set_title('Average Nonlinear Iteration Count per Timestep',fontsize=16)
ax2.set_xlabel('dt',fontsize=16)
ax2.set_ylim(0,avgPlotMax+1)
ax2.legend(loc='best',fontsize=16)
ax2.tick_params(axis='both', which='major', labelsize=14)


pyplot.show()
