import matplotlib.pyplot as pyplot
import numpy as np
import matplotlib

matplotlib.rcParams.update({'font.size': 24})

fig, axList = pyplot.subplots(1,3,figsize=(30,10))
axQl = axList[0]
axQv = axList[1]
axAl = axList[2]

f = 0.3
qsat = 1.0
ymin = -0.01
ymax = 1.01
# background
for ax in axList:
  ax.fill([0.0, f, f, 0.0], [ymin, ymin, ymax, ymax], color=(0.0,0.0,0.0,0.2))
  ax.fill([f, 1.0, 1.0, f], [ymin, ymin, ymax, ymax], color=(0.0,1.0,1.0,0.3))

# ql
axQl.plot([f, 1.0], [0.0, 0.0], '-k', linewidth=3)
ql = f/2.0
x = np.linspace(0, f, 100)
a0 = ql/f
a1 = -2.0*ql/f**2
a2 = 3.0*ql/f**3
a4 = 5.0*ql/f**5
axQl.plot(x, a0*np.ones(np.shape(x)), '-m', label='constant', linewidth=3)
axQl.plot([f, f], [0, a0], ':m', linewidth=7)
axQl.plot(x, a1*(x-f), '--r', label='linear', linewidth=3)
axQl.plot(x, a2*(x-f)**2, '-.g', label='quadratic', linewidth=3)
axQl.plot(x, a4*(x-f)**4, '-b', label='quartic', linewidth=3)
axQl.set_ylabel('liquid concentration', fontsize='x-large')

# qv
axQv.plot([0, f], [qsat, qsat], '-k', linewidth=3)
qv = (1.0-f)/2.0
x = np.linspace(f, 1, 100)
a0 = (qv - (1.0-f)*qsat)/(1.0-f)
a1 = 2.0*(qv - (1.0-f)*qsat)/(1.0-f)**2
a2 = 3.0*(qv - (1.0-f)*qsat)/(1.0-f)**3
a4 = 5.0*(qv - (1.0-f)*qsat)/(1.0-f)**5
axQv.plot(x, qsat + a0*np.ones(np.shape(x)), '-m', label='constant', linewidth=3)
axQv.plot([f,f], [qsat + a0, qsat], ':m', linewidth=7)
axQv.plot(x, qsat + a1*(x-f), '--r', label='linear', linewidth=3)
axQv.plot(x, qsat + a2*(x-f)**2, '-.g', label='quadratic', linewidth=3)
axQv.plot(x, qsat + a4*(x-f)**4, '-b', label='quartic', linewidth=3)
axQv.set_ylabel('vapor concentration', fontsize='x-large')

# Al
Al = 1.0/2.0
x = np.linspace(0, 1, 100)
a0 = Al
a1 = -2.0*Al
a2 = 3.0*Al
a4 = 5.0*Al
axAl.plot(x, a0*np.ones(np.shape(x)), '-m', label='constant', linewidth=3)
axAl.plot(x, a1*(x-1.0), '--r', label='linear', linewidth=3)
axAl.plot(x, a2*(x-1.0)**2, '-.g', label='quadratic', linewidth=3)
axAl.plot(x, a4*(x-1.0)**4, '-b', label='quartic', linewidth=3)
axAl.set_ylabel('liquid concentration tendency', fontsize='x-large')


for ax in axList:
  ax.set_xlim(0,1)
  ax.set_ylim(ymin, ymax)
  ax.legend(loc='best', fontsize='large')
  ax.set_xlabel('x', fontsize='x-large')

fig.tight_layout()
fig.savefig('figure1-profiles.pdf')

pyplot.show()