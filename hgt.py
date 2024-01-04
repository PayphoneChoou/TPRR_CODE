import xarray as xr
import numpy as np
import netCDF4 as nc
import time as tm


dfile = '/raid63/WRF_reinit_GLDASsnow_SKINsst/Output/2019/2019010112/wrfout_d01_2019-01-02_12:00:00'
data_reinit_sn = xr.open_dataset(dfile)
lon_reinit_sn,lat_reinit_sn = data_reinit_sn.XLONG[0].values,data_reinit_sn.XLAT[0].values


nc_file_path = '/raid63/pfzhou/DATA/WRF_Reinit_SN/nc/static/WRF_reinit_9km_static_hgt.nc'
root = nc.Dataset(nc_file_path,mode = 'w',format='NETCDF4')
# 为nc文件创建维度
lon = root.createDimension('lon',530)
lat = root.createDimension('lat',360)

# 为nc文件创建变量
longitude = root.createVariable('longitude',np.float32,('lat','lon'))
latitude = root.createVariable('latitude',np.float32,('lat','lon'))
hgt = root.createVariable('hgt',np.float32,('lat','lon'))
#SWDOWN,GLW
# 对nc文件进行描述
root.description = 'Nanjing University,Science of Atmosphere'
root.history = 'Created ' + tm.ctime(tm.time())
root.source = 'netCDF4 python '
latitude.units = 'degree_nonth'
longitude.units = 'degree_east'
hgt.description = 'terrain height'

longitude[:,:] = lon_reinit_sn
latitude[:,:] = lat_reinit_sn
hgt[:,:] = data_reinit_sn.HGT.values

root.close()

