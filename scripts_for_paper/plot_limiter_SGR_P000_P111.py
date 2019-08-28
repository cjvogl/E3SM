import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import os.path

matplotlib.rcParams.update({'font.size': 24})

fileDict = {'SGR-P[0,0,0]': './data/SGR_P000_RKZ_qlneg_lm4_855579.txt',
            'SGR-P[1,1,1]': './data/SGR_P111_RKZ_qlneg_lm4_855579.txt'}
linespecDict = {'SGR-P[0,0,0]': '--g',
                'SGR-P[1,1,1]': '-b'}

fig, ax = pyplot.subplots(figsize=(10,10))
for sgr in fileDict.keys():
  tmp = np.loadtxt(fileDict[sgr])
  t = tmp[0,:]
  val = tmp[1,:]
  ax.plot(t[2:], val[2:], linespecDict[sgr], label=sgr, linewidth=3)

ax.set_ylabel('clipped condensation (kg/kg)', fontsize='x-large')
ax.set_xlabel('time (s)', fontsize='x-large')
ax.legend(loc='best', fontsize='large')
fig.tight_layout()

fig.savefig('clipping_limiter.pdf')
pyplot.show()
