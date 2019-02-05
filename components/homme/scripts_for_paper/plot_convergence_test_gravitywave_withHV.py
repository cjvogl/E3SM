import matplotlib.pyplot as pyplot
import matplotlib
from dictionaries import methodDict, lineStyleDict, colorDict
import numpy as np

matplotlib.rcParams.update({'font.size':22})

methodList = ['KGU35', 'ARS232', 'DBM453', 'ARS222', 'ARS233', 'ARS343', \
              'ARS443', 'ARK324', 'ARK436', 'SSP3333b', 'SSP3333c']
epsilon_times_1000 = 2.220446049250313e-13
roundoff = epsilon_times_1000*768.0

# iterate through methods and plot convergence test
f, ax = pyplot.subplots(figsize=(10,10))
for method in methodList:
  tmp = np.loadtxt('./data/' + method + '_convergence_data_withHV.txt')
  dtPlot = tmp[0,:]
  errorPlot = tmp[1,:]
  # compute order as slope between the smallest 2 timesteps where temporal error
  # is dominating instead of roundoff error
  ind = 0
  while (errorPlot[ind] < roundoff):
    ind = ind + 1
  order = np.log(errorPlot[ind+1]/errorPlot[ind])/np.log(dtPlot[ind+1]/dtPlot[ind])
  # plot
  ax.loglog(dtPlot, errorPlot, lineStyleDict[method], \
            label='%s (%3.2f)' % (method,order), \
            color=colorDict[method], linewidth=3, markersize=12)
  ax.loglog(dtPlot, roundoff*np.ones(np.shape(dtPlot)), '--k', label=None, linewidth=2)

# finalize plots and save without legend but show with legend
ax.set_ylabel('maximum relative error', fontsize='xx-large')
ax.set_xlabel('dt (s)', fontsize='xx-large')
ax.set_xlim(6e-4,7e0)
ax.set_ylim(5e-11,1e-6)
#ax.axis('equal')
f.tight_layout()
f.savefig('convergence_test_gravitywave_noHV.pdf')
ax.legend(loc='best')
pyplot.show()
