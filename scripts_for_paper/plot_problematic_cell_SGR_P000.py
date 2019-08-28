import matplotlib
import matplotlib.pyplot as pyplot
import numpy as np
import os.path

matplotlib.rcParams.update({'font.size': 24})

directory = 'RKZ_SGR_extrapf_qv0_ql0_Al0_lmt4_adjIC_ne30_ne30_FC5AQUAP_intel_quartz/DT0004_01_dycore_mac/'
suffix = '_855579.txt'

txtFile = './data/SGR_P000_RKZ_term_A_855579.txt'
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
term_A = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_term_B_855579.txt'
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
term_B = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_term_C_855579.txt'
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
term_C = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_qme_855579.txt'
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
Q = tmp[1,:]

txtFile = './data/SGR_P000_RKZ_fac_855579.txt'
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
f = tmp[1,:]

txtFile = './data/SGR_P000_CLDLIQ_855579.txt'
tmp = np.loadtxt(txtFile)
t = tmp[0,:]
ql = tmp[1,:]

fig, ax = pyplot.subplots(figsize=(11,10))
ax.plot(t[2:], ql[2:], '-k', label='cloud liquid', linewidth=3)
ax.plot(t[2:], term_A[2:] + term_B[2:] + term_C[2:], '--g', label='original Q', linewidth=3)
ax.plot(t[2:], Q[2:], ':r', label='clipped Q', linewidth=3)
ax2 = ax.twinx()
ax2.plot(t[2:], f[2:], '-.b', linewidth=3)
ax.set_xlabel('time (s)', fontsize='x-large')
ax.set_ylabel('concentration (kg/kg)', fontsize='x-large')
ax2.set_ylabel('cloud fraction', color='b', fontsize='x-large')
ax2.tick_params('y', colors='b')
ax.legend(loc='upper center', fontsize='large')
fig.tight_layout()

fig.savefig('problematic_cell_P000.pdf')
pyplot.show()
