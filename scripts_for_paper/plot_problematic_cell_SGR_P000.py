import matplotlib
import matplotlib.pyplot as pyplot
from matplotlib.lines import Line2D
import numpy as np
import os.path

matplotlib.rcParams.update({'font.size': 24})

dt = 4.0
directory = 'RKZ_SGR_extrapf_qv0_ql0_Al0_lmt4_adjIC_ne30_ne30_FC5AQUAP_intel_quartz/DT0004_01_dycore_mac/'
suffix='_1069852.txt'

txtFile = './data/SGR_P000_RKZ_term_A' + suffix
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
term_A = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_term_B' + suffix
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
term_B = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_term_C' + suffix
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
term_C = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_qme' + suffix
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
Q = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_fac' + suffix
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
f = tmp[1,:]

txtFile = './data/SGR_P000_CLDLIQ' + suffix
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
ql = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_Al' + suffix
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
Al = tmp[1,:]

fig, ax = pyplot.subplots(1, 2, sharex=True, sharey=True, figsize=(20,10))

ax[0].plot(t[2:], dt*(term_A[2:] + term_B[2:] + term_C[2:]), '-g', label='dt*(A+B+C)', linewidth=3)
ax[0].plot(t[2:], dt*term_C[2:], '-r', label='dt*C', linewidth=3)
ax2 = ax[0].twinx()
ax2.plot(t[2:], f[2:], '-b', linewidth=3)
ax2.tick_params('y', colors='w')
ax[0].set_xlabel('time (s)', fontsize='x-large')
ax[0].set_ylabel('concentration (kg/kg)', fontsize='x-large')
ax[0].set_ylim(-4.0e-8, 3.0e-8)
handles, labels = ax[0].get_legend_handles_labels()
handles.append(Line2D([], [], color='blue', linewidth=3, linestyle='-'))
labels.append('cloud frac.')
ax[0].legend(handles, labels, loc='lower right', fontsize='large')

ax[1].plot(t[2:], ql[2:], '-k', label='cloud liquid', linewidth=3)
ax[1].plot(t[2:], dt*Al[2:], '--m', label='dt*(liq. tend.)', linewidth=3)
ax2 = ax[1].twinx()
ax2.plot(t[2:], f[2:], '-b', linewidth=3)
ax2.tick_params('y', colors='b')
ax2.set_ylabel('fraction of grid cell', color='b', fontsize='x-large')
ax[1].set_xlabel('time (s)', fontsize='x-large')
handles, labels = ax[1].get_legend_handles_labels()
handles.append(Line2D([], [], color='blue', linewidth=3, linestyle='-'))
labels.append('cloud frac.')
ax[1].legend(handles, labels, loc='lower right', fontsize='large')

fig.tight_layout()

fig.savefig('problematic_cell_P000.pdf')
pyplot.show()
