import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import pdb
import os.path
#
path = '/Volumes/Zhuomin8T/MAPP/EN4'
# year range: 1954-2024
yr1, yr2 = [1954, 2024]
ny = yr2 - yr1 + 1
# lat range:
a1, a2 = [118, 149]
# lon range:
b1, b2 = [279, 320]
# depth range around 200 m
c1, c2 = [15, 17]
#
T200 = np.ma.zeros([ny,12,2,a2-a1,b2-b1]);T200.mask = True
T207 = np.ma.zeros([ny,12,a2-a1,b2-b1]);T207.mask = True
T185 = T207.copy()
for year in range(yr1,yr2+1):
    for mon in range(12):
        filename = '%s/EN.4.2.2.f.analysis.g10.%d%02d.nc' % (path,year,mon+1)
        if os.path.exists(filename):
            nc = netCDF4.Dataset(filename,'r')
            lat = nc.variables['lat'][a1:a2]
            lon = nc.variables['lon'][b1:b2] - 360.
            depth = nc.variables['depth'][c1:c2]
            T200[year-yr1,mon] = nc.variables['temperature'][0,c1:c2,a1:a2,b1:b2] - 273.15# change unit to degC
            T185[year-yr1,mon] = nc.variables['temperature'][0,c1,a1:a2,b1:b2] - 273.15# change unit to degC
            T207[year-yr1,mon] = nc.variables['temperature'][0,c1+1,a1:a2,b1:b2] - 273.15# change unit to degC
            nc.close()
#
print(depth)
print(lat.min(),lat.max())
print(lon.min(),lon.max())
nlat = len(lat)
nlon = len(lon)
dlon,dlat = np.meshgrid(lon,lat)
T200 = np.reshape(T200,[ny*12,2,nlat,nlon])
T185 = np.reshape(T185,[ny*12,nlat,nlon])
T207 = np.reshape(T207,[ny*12,nlat,nlon])
#T200 = np.mean(T200,axis=0) # 2,31,41
# 9 points
# interpolation to 200m depth
T200_m = T200[:,0] + (200-depth[0])*(T200[:,1]-T200[:,0])/(depth[1]-depth[0]) # ny*12,31,41
## save to file
output = '/Users/zhuominchen/Documents/Postdoc/Prediction/MAPP/CMEMS/DIAGS'
file1 = '%s/EN4_T200_monthly_%d-%d.nc' % (output,yr1,yr2)
root = netCDF4.Dataset(file1,'w',format='NETCDF4')
#root.createDimension('num',len(po_lon))
root.createDimension('nmon',ny*12)
root.createDimension('nlat',nlat)
root.createDimension('nlon',nlon)
##
tT185 = root.createVariable('T185','f4',('nmon','nlat','nlon'))
tT200 = root.createVariable('T200','f4',('nmon','nlat','nlon'))
tT207 = root.createVariable('T207','f4',('nmon','nlat','nlon'))
tlat = root.createVariable('lat','f4',('nlat',))
tlon = root.createVariable('lon','f4',('nlon',))
tT200[:,:,:] = T200_m[:,:,:]
tT185[:,:,:] = T185[:,:,:]
tT207[:,:,:] = T207[:,:,:]
tlat[:] = lat[:]
tlon[:] = lon[:]
root.close()
