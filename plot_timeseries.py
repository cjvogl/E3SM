from netCDF4 import Dataset
from matplotlib.pyplot import colorbar, subplots, show
from mpl_toolkits.basemap import Basemap
import numpy as np

timeList = np.array([0, 1800, 3600, 5400, 7200])
#config = 'RKZ_P3_lmt4_adjIC'
config = 'RKZ_P3_lmt4'
#timeList = [0, 1800, 3600]
#config = 'RKZ_A1_B1_C2_ql19_lmt4_adjIC'
f, ax = subplots(6, len(timeList), sharex=True, sharey=True, figsize=(16,10))
f2, ax2 = subplots(3, len(timeList), figsize=(16,10))
f3, ax3 = subplots(3, len(timeList), figsize=(16,10))

# create objects for extrema
data = Dataset('./%s_ne30_ne30_FC5AQUAP_intel_cori-knl/DT0001_01_dycore_mac/DT0001_01_dycore_mac.cam.h0.0001-01-01-%05d.nc' % (config, timeList[-1]))
ql = data['CLDLIQ'][:]
qv = data['Q'][:]
T = data['T'][:]
p0 = data['P0'][:]
ps = data['PS'][:]
ps = ps[0,:]
a = data['hyai'][:]
b = data['hybi'][:]
dp = np.outer(a[1:] - a[0:-1], p0*np.ones(np.shape(ps))) + np.outer(b[1:] - b[0:-1], ps)
qlTotal = np.sum(ql[0,:,:]*dp, axis=0)
qvTotal = np.sum(qv[0,:,:]*dp, axis=0)
Ttotal = np.sum(T[0,:,:]*dp, axis=0)
indQl = np.argmax(qlTotal)

# create objects to store previous solution
data = Dataset('./%s_ne30_ne30_FC5AQUAP_intel_cori-knl/DT0001_01_dycore_mac/DT0001_01_dycore_mac.cam.h0.0001-01-01-00000.nc' % config)
ql = data['CLDLIQ'][:]
qv = data['Q'][:]
T = data['T'][:]
p0 = data['P0'][:]
ps = data['PS'][:]
ps = ps[0,:]
a = data['hyai'][:]
b = data['hybi'][:]
dp = np.outer(a[1:] - a[0:-1], p0*np.ones(np.shape(ps))) + np.outer(b[1:] - b[0:-1], ps)
qlPrev = ql[0,:,indQl]
qlTotalPrev = np.sum(ql[0,:,:]*dp, axis=0)
qvPrev = qv[0,:,indQl]
qvTotalPrev = np.sum(qv[0,:,:]*dp, axis=0)
Tprev = T[0,:,indQl]
TtotalPrev = np.sum(T[0,:,:]*dp, axis=0)


# plot scaling
minQl = 0.0
maxQl = 100.0
dQl = 20.0
minT = 221.0
maxT = 260.0
dT = 0.2
minQv = 4.0
maxQv = 700.0
dQv = 20.0

for i, time in enumerate(timeList):
  data = Dataset('./%s_ne30_ne30_FC5AQUAP_intel_cori-knl/DT0001_01_dycore_mac/DT0001_01_dycore_mac.cam.h0.0001-01-01-%05d.nc' % (config, time))
  ql = data['CLDLIQ'][:]
  qv = data['Q'][:]
  T = data['T'][:]
  lon = data['lon'][:]
  lat = data['lat'][:]
  p0 = data['P0'][:]
  ps = data['PS'][:]
  ps = ps[0,:]
  a = data['hyai'][:]
  b = data['hybi'][:]
  dp = np.outer(a[1:] - a[0:-1], p0*np.ones(np.shape(ps))) + np.outer(b[1:] - b[0:-1], ps)
  a = data['hyam'][:]
  b = data['hybm'][:]
  p = np.outer(a, p0*np.ones(np.shape(ps))) + np.outer(b, ps)
  qlTotal = np.sum(ql[0,:,:]*dp, axis=0)
  qvTotal = np.sum(qv[0,:,:]*dp, axis=0)
  Ttotal = np.sum(T[0,:,:]*dp, axis=0)
  Z = np.sum(dp,axis=0)
  
  #### Figure 3 ####
  x = np.outer(np.ones(30), [0.0,1.0])
  z = np.outer((ps[indQl]-p[:,indQl])/9.8, [1.0,1.0])/1000.0
  qlPlot = np.outer(ql[0,:,indQl], [1.0,1.0])
  qlDiffPlot = np.outer(ql[0,:,indQl]-qlPrev, [1.0,1.0])
  tmp = ax3[0,i].pcolor(x, z, qlDiffPlot, cmap='PiYG', vmin=-0.001, vmax=0.001)
  colorbar(tmp, ax=ax3[0,i])
  ax3[0,i].set_xticks([])
  ax3[0,i].set_title('vert. prof. @ %ds' % time)
  ax3[0,i].set_ylabel('approx. elev. (km)')
  qlPrev = ql[0,:,indQl]
  qvPlot = np.outer(qv[0,:,indQl], [1.0,1.0])
  qvDiffPlot = np.outer(qv[0,:,indQl]-qvPrev, [1.0,1.0])
  tmp = ax3[1,i].pcolor(x, z, qvDiffPlot, cmap='PiYG', vmin=-0.001, vmax=0.001)
  colorbar(tmp, ax=ax3[1,i])
  ax2[1,i].set_xticks([])
  ax2[1,i].set_ylabel('approx. elev. (km)')
  qvPrev = qv[0,:,indQl]
  TdiffPlot = np.outer(T[0,:,indQl]-Tprev, [1.0,1.0])
  tmp = ax3[2,i].pcolor(x, z, TdiffPlot, cmap='bwr', vmin=-5.0, vmax=5.0)
  colorbar(tmp, ax=ax3[2,i])
  ax2[2,i].set_xticks([])
  ax2[2,i].set_ylabel('approx. elev. (km)')
  Tprev = T[0,:,indQl]

  #### Figure 2 ####
  ax2[0,0].plot(time, qlTotal[indQl],'ko')
  ax2[1,0].plot(time, qvTotal[indQl],'ko')
  ax2[2,0].plot(time, Ttotal[indQl]/Z[indQl],'ko')
  if (i > 0):
    x = np.outer(np.ones(30), [0.0,1.0])
    z = np.outer((ps[indQl]-p[:,indQl])/9.8, [1.0,1.0])/1000.0
    tmp = ax2[0,i].pcolor(x, z, qlPlot, cmap='Greens', vmin=0.0, vmax=0.0035)
    qlPlot = np.outer(ql[0,:,indQl], [1.0,1.0])
    colorbar(tmp, ax=ax2[0,i])
    ax2[0,i].set_xticks([])
    ax2[0,i].set_title('vert. diff. prof. @ %ds' % time)
    ax2[0,i].set_ylabel('approx. elev. (km)')
    qvPlot = np.outer(qv[0,:,indQl], [1.0,1.0])
    tmp = ax2[1,i].pcolor(x, z, qvPlot, cmap='Greens', vmin=0.0, vmax=0.013)
    colorbar(tmp, ax=ax2[1,i])
    ax2[1,i].set_xticks([])
    ax2[1,i].set_ylabel('approx. elev. (km)')
    Tplot = np.outer(T[0,:,indQl], [1.0,1.0])
    tmp = ax2[2,i].pcolor(x, z, Tplot, cmap='jet', vmin=210, vmax=290)
    colorbar(tmp, ax=ax2[2,i])
    ax2[2,i].set_xticks([])
    ax2[2,i].set_ylabel('approx. elev. (km)')

  #### Figure 1 ####

  # plot total cloud liquid
  bmap = Basemap(lon_0=180)
  x,y = bmap(lon,lat)
  bmap.ax = ax[0,i]
  sQl = bmap.scatter(x, y, c=qlTotal, cmap='Greens', vmin=minQl, vmax=maxQl)
  print(np.amin(qlTotal),np.amax(qlTotal))
  ax[0,i].set_title('%ds' % time)
  bmap.drawcoastlines()
  # plot total water vapor
  bmap = Basemap(lon_0=180)
  x,y = bmap(lon,lat)
  bmap.ax = ax[1,i]
  sQv = bmap.scatter(x, y, c=qvTotal, cmap='Greens', vmin=minQv, vmax=maxQv)
  print(np.amin(qvTotal),np.amax(qvTotal))
  bmap.drawcoastlines()
  # plot average temperature
  bmap = Basemap(lon_0=180)
  x,y = bmap(lon,lat)
  bmap.ax = ax[2,i]
  sT = bmap.scatter(x, y, c=Ttotal/Z, cmap='jet', vmin=minT, vmax=maxT)
  print(np.amin(Ttotal/Z),np.amax(Ttotal/Z))
  bmap.drawcoastlines()
  # plot difference in total cloud liquid
  bmap = Basemap(lon_0=180)
  x,y = bmap(lon,lat)
  bmap.ax = ax[3,i]
  sQlDiff = bmap.scatter(x, y, c=qlTotal-qlTotalPrev, cmap='PiYG', vmin=-dQl, vmax=dQl)
  bmap.drawcoastlines()
  qlTotalPrev = qlTotal
  # plot average difference in total water vapor
  bmap = Basemap(lon_0=180)
  x,y = bmap(lon,lat)
  bmap.ax = ax[4,i]
  sQvdiff = bmap.scatter(x, y, c=(qvTotal-qvTotalPrev), cmap='PiYG', vmin=-dQv, vmax=dQv)
  bmap.drawcoastlines()
  qvTotalPrev = qvTotal
  # plot average difference in temperature
  bmap = Basemap(lon_0=180)
  x,y = bmap(lon,lat)
  bmap.ax = ax[5,i]
  sTdiff = bmap.scatter(x, y, c=(Ttotal-TtotalPrev)/Z, cmap='bwr', vmin=-dT, vmax=dT)
  bmap.drawcoastlines()
  TtotalPrev = Ttotal

ax[0,0].set_ylabel('Total CLDLIQ')
ax[1,0].set_ylabel('Total Q')
ax[2,0].set_ylabel('Avg Temp')
ax[3,0].set_ylabel('Total CLDLIQ Diff')
ax[4,0].set_ylabel('Total Q Diff')
ax[5,0].set_ylabel('Avg Temp Diff')


ax2[0,0].set_ylabel('Total CLDLIQ')
ax2[1,0].set_ylabel('Total Q')
ax2[2,0].set_ylabel('Avg Temp')
ax2[2,0].set_xlabel('time (s)')


f.tight_layout()
f2.tight_layout()

# add colorbars
width = 0.01
right = 0.92
f.subplots_adjust(right=right)
bbox = ax[0,0].get_position()
cax = f.add_axes([right+0.01, bbox.p0[1], width, bbox.height])
colorbar(sQl, cax=cax)
bbox = ax[1,0].get_position()
cax = f.add_axes([right+0.01, bbox.p0[1], width, bbox.height])
colorbar(sT, cax=cax)
bbox = ax[2,0].get_position()
cax = f.add_axes([right+0.01, bbox.p0[1], width, bbox.height])
colorbar(sQv, cax=cax)
bbox = ax[3,0].get_position()
cax = f.add_axes([right+0.01, bbox.p0[1], width, bbox.height])
colorbar(sQlDiff, cax=cax)
bbox = ax[4,0].get_position()
cax = f.add_axes([right+0.01, bbox.p0[1], width, bbox.height])
colorbar(sQvdiff, cax=cax)
bbox = ax[5,0].get_position()
cax = f.add_axes([right+0.01, bbox.p0[1], width, bbox.height])
colorbar(sTdiff, cax=cax)
show()

 
  
