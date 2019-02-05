import matplotlib.pyplot as pyplot
import matplotlib
from dictionaries import methodDict
import numpy as np

matplotlib.rcParams.update({'font.size':22})

method = 'ARS343'

f, ax = pyplot.subplots(figsize=(10,10))
lineStyles = ['x', '-']
tsteps = [10, 20]
colors = ['tab:blue', 'tab:green', 'tab:red', 'k']

for l,tstep in enumerate(tsteps):
  fileNames = ['./data/tsteptype%d_tstep%d_restart64800_ndays1_nu0.0_rsplit0_energy_error_data.txt' % (methodDict[method], tstep),
               './data/tsteptype%d_tstep%d_restart64800_ndays1_nu1e15_rsplit0_energy_error_data.txt' % (methodDict[method], tstep),
               './data/tsteptype%d_tstep%d_restart64800_ndays1_nu0.0_rsplit1_energy_error_data.txt' % (methodDict[method], tstep),
               './data/tsteptype%d_tstep%d_restart64800_ndays1_nu1e15_rsplit1_energy_error_data.txt' % (methodDict[method], tstep)]
  for m,fileName in enumerate(fileNames):
    tmp = np.loadtxt(fileName)
    timeList = tmp[0,:]
    errorList = tmp[1,:]
    ax.plot(timeList, errorList, lineStyles[l], color=colors[m],
            linewidth=3, markersize=12)
ax.set_xlabel('time (h)', fontsize='xx-large')
ax.set_ylabel('relative conservation error', fontsize='xx-large')
ax.set_xlim(0,24)
ax.set_ylim(-1.5e-7,1.5e-7)

f.tight_layout()
f.savefig('energy_error_baroclinicinstability_ars343.pdf')
pyplot.show()
