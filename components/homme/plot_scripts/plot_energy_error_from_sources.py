import glob
from matplotlib.pyplot import plot, show, subplots
import matplotlib
from dictionaries import methodDict, lineStyleDict, colorDict

matplotlib.rcParams.update({'font.size':16})

f, ax = subplots(2,2, sharex=True, sharey=True, figsize=(10,10))
f2, ax2 = subplots(figsize=(10,10))
suffixList = ('nu0.0_rsplit0', \
              'nu1e15_rsplit0', \
              'nu0.0_rsplit1', \
              'nu1e15_rsplit1')

for j,suffix in enumerate(suffixList):
  for method in methodDict.keys():
    globstr = 'tsteptype%d_tstep10_restart64800_ndays1_%s.out' % (methodDict[method],suffix)
    for fileName in glob.glob(globstr):
      print('Reading in from ' + fileName)
      words = fileName.split('_')
      dt = words[1].replace('tstep','')
  #    dt = dt.replace('.out','')
      fileObj = open(fileName)
      lines = list(fileObj)
      fileObj.close()
      time = 0.0
      timeList = [time]
      errorList = [0.0,]
      for line in lines:
        if ('time=' in line):
          time = float(line.split()[3])
          if ('[s]' in line):
            time /= 3600.0
        if ('(E-E0)/E0' in line):
          timeList.append(time)
          errorList.append(float(line.split()[1]))
      ax[j/2,j%2].plot(timeList, errorList, lineStyleDict[method], label=method, color=colorDict[method])
      if (j == 0):
        ax2.plot(timeList, errorList, lineStyleDict[method], label=method, color=colorDict[method])
  ax[j/2,j%2].set_title(suffix)
  if (j/2 == 1):
    ax[1,j%2].set_xlabel('time (h)')
  if (j%2 == 0):
    ax[j/2,0].set_ylabel('global energy cons error')
  if (j == 0):
    ax2.set_xlabel('time (h)')
    ax2.set_ylabel('global energy cons error')
    #ax2.set_ylim(-6.0e-13,7.0e-13)
    ax2.legend(loc='best')

f.tight_layout()
f2.tight_layout()
show()
