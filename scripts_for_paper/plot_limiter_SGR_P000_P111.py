import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import os.path

matplotlib.rcParams.update({'font.size': 24})

directoryDict = {'SGR-P[0,0,0]': 'RKZ_SGR_extrapf_qv0_ql0_Al0_lmt4_adjIC_ne30_ne30_FC5AQUAP_intel_quartz/DT0004_01_dycore_mac/',
                 'SGR-P[1,1,1]': 'RKZ_SGR_extrapf_qv1_ql1_Al1_lmt4_adjIC_ne30_ne30_FC5AQUAP_intel_quartz/DT0004_01_dycore_mac/'}
linespecDict = {'SGR-P[0,0,0]': '--g',
                'SGR-P[1,1,1]': '-b'}
suffix = '_855579.txt'

fig, ax = pyplot.subplots(figsize=(10,10))
for sgr in directoryDict.keys():
  txtFile = directoryDict[sgr] + 'RKZ_qlneg_lm4' + suffix
  tmp = np.loadtxt(txtFile)
  t = tmp[0,:]
  val = tmp[1,:]
  ax.plot(t[2:], val[2:], linespecDict[sgr], label=sgr, linewidth=3)

ax.set_ylabel('clipped condensation (kg/kg)', fontsize='x-large')
ax.set_xlabel('time (s)', fontsize='x-large')
ax.legend(loc='best', fontsize='large')
fig.tight_layout()

fig.savefig('clipping_limiter.pdf')
pyplot.show()
