import glob
import matplotlib.pyplot as pyplot
from netCDF4 import Dataset
import numpy as np

noHV = False

methodDict = {#'U35-ref': 5,
                  #'U35': 11,
                  #'ARS232-ref': 7,
                  #'ARS232': 12,
                  'DBM453': 13}
                  #'ARS222': 14,
                  #'ARS233': 15,
                  #'ARS343': 16,
                  #'ARS443': 17,
                  #'ARK324': 18}
                  #'ARK436': 19}

# Iterate through methods and plot nonlinear iteration count
f, ax1 = pyplot.subplots(figsize=(10,10))
for m,method in enumerate(methodDict.keys()):
  dtList = []
  countDict = {}
  print method
  # Load timestep and solution data
  for fileName in glob.glob('tsteptype%d_tstep*.out' % methodDict[method]):
    if ((noHV and "nu0.0" in fileName) or (not noHV and "nu0.0" not in fileName)):
      words = fileName.split('_')
      dt = words[1].replace('tstep','')
      dtList.append(float(dt))
      # obtain nonlinear iteration count
      f = open(fileName)
      lines = list(f)
      f.close
      for line in reversed(lines):
        if (line.startswith('  Max num nonlin iters')):
          words = line.split()
          countDict[dt] = int(words[5])
          break
  # Sort values and add to plot
  dtPlot = np.sort(dtList)
  countPlot = np.empty(len(dtPlot))
  for j,dt in enumerate(dtPlot):
    countPlot[j] = countDict[str(dt)]
  if (m < 7):
    lineStyle = '-o'
  else:
    lineStyle = '-d'
  ax1.plot(dtPlot,countPlot,lineStyle,label='%s' % method)

#orderList = ('first', 'second', 'third')
#for j, order in enumerate(orderList):
#    coeff = L2Plot[0]/dtPlot[0]**(j+1)
#    ax1.loglog(dtPlot,coeff*dtPlot**(j+1),'--',label=order)
#    coeff = LIPlot[0]/dtPlot[0]**(j+1)
#    ax2.loglog(dtPlot,coeff*dtPlot**(j+1),'--',label=order)

ax1.set_title('Average Nonlinear Iteration Count per Timestep')
ax1.set_xlabel('dt')
ax1.legend(loc='best')

pyplot.show()
