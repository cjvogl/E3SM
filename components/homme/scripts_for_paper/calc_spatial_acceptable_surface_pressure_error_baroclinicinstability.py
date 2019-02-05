from netCDF4 import Dataset
from mpl_toolkits.basemap import Basemap
import numpy as np
import sys

if (len(sys.argv) < 2):
  print "Usage: python calc_spatial_acceptable_surface_pressure_error.py <base_output_directory>"
  exit()
outpath = sys.argv[1]

# read in reference solution
fileRef = outpath + '/output_tsteptype5_tstep10_hydrostatic/dcmip2012_test41.nc'
print 'Reading reference solution from ' + fileRef
data = Dataset(fileRef)
qRef = data['ps'][:]
lon = data['lon'][:]
lat = data['lat'][:]

# read in acceptable solution
fileRef = outpath + '/output_tsteptype5_tstep300_hydrostatic/dcmip2012_test41.nc'
print 'Reading acceptable solution from ' + fileRef
data = Dataset(fileRef)
t = data['time'][:]
q = data['ps'][:]
numElements = np.shape(q)[1]*np.shape(q)[2]
errorTol = np.sqrt( np.sum((q-qRef)**2, axis=(1,2))/numElements )

bmap = Basemap(lon_0=180)
x,y = bmap(lon, lat)
X = np.outer(np.ones(np.shape(y)), x)
Y = np.outer(y, np.ones(np.shape(x)))

output = np.vstack((X, Y, (qRef[-1,:,:] - qRef[0,:,:])/1e3))
filename = './data/surface_pressure_reference_data.txt'
np.savetxt(filename, output)

output = np.vstack((X, Y, (q[-1,:,:] - qRef[-1,:,:])/1e3))
filename = './data/surface_pressure_tolerance_spatial_data.txt'
np.savetxt(filename, output)

output = np.vstack((t, errorTol/1e3))
filename = './data/surface_pressure_tolerance_rms_data.txt'
np.savetxt(filename, output)
