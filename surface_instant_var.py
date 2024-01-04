import os
import xarray as xr
import numpy as np
import pandas as pd
import datetime
import netCDF4 as nc 
import time as tm
import sys


st_year,ed_year,more,dpath1,dpath2,var,unit,var_file = int(sys.argv[1]),int(sys.argv[2]),int(sys.argv[3]),str(sys.argv[4]),str(sys.argv[5]),str(sys.argv[6]),str(sys.argv[7]),str(sys.argv[8])
print('start year = ',st_year,'\n','end_year = ',ed_year,'\n','var = ',var,'unit = ',unit,'var_file = ',var_file)
dpath3 = '/raid63/pfzhou/DATA/'
for year in range(st_year,ed_year,1):
    dfile = '/raid63/WRF_reinit_GLDASsnow_SKINsst/Output/2019/2019010112/wrfout_d01_2019-01-01_12:00:00'
    data_reinit_sn = xr.open_dataset(dfile)
    lon_reinit_sn,lat_reinit_sn = data_reinit_sn.XLONG[0].values,data_reinit_sn.XLAT[0].values
    
    #%%
    dtimes = pd.date_range(start='{}-12-31 12:00:00'.format(year-1),end='{}-12-30 12:00:00'.format(year),freq='24H').strftime('%Y%m%d%H')
    
    dlists = []

    for i in range(more):
        dlists.append(dpath1 + '/{}/{}'.format(year-1,dtimes[i]))    
    for i in range(more,len(dtimes)):
        dlists.append(dpath2 + '/{}/{}'.format(year,dtimes[i]))
     
    nc_file_path = dpath3 + 'WRF_reinit_9km_h_sfc_{}_{}.nc'.format(var_file,year)
    root = nc.Dataset(nc_file_path,mode = 'w',format='NETCDF4')
    # 为nc文件创建维度
    lon = root.createDimension('lon',530)
    lat = root.createDimension('lat',360)
    time = root.createDimension('time',None)
    
    # 为nc文件创建变量
    times = root.createVariable('times',int,('time'))
    longitude = root.createVariable('longitude',np.float32,('lat','lon'))
    latitude = root.createVariable('latitude',np.float32,('lat','lon'))
    
    VAR = root.createVariable(var,np.float32,('time','lat','lon'))
    #SWDOWN,GLW
    # 对nc文件进行描述
    root.description = 'Nanjing University,Science of Atmosphere'
    root.history = 'Created ' + tm.ctime(tm.time())
    root.source = 'netCDF4 python '
    latitude.units = 'degree_nonth'
    longitude.units = 'degree_east'
    VAR.units = unit
    
    longitude[:,:] = lon_reinit_sn
    latitude[:,:] = lat_reinit_sn
    
    count = 0
    for i,dlist in enumerate(dlists):
        if os.path.isdir(dlist):
            print(i,dlist,dtimes[i])
            start_time = (datetime.datetime.strptime(dtimes[i],'%Y%m%d%H')+datetime.timedelta(hours=12)).strftime('%Y %m %d %H')
            end_time = (datetime.datetime.strptime(dtimes[i],'%Y%m%d%H')+datetime.timedelta(hours=36)).strftime('%Y %m %d %H')
            dhours = pd.date_range(start=start_time,end=end_time,freq='H').strftime('%Y-%m-%d_%H:%M:%S')
            hours = pd.date_range(start=start_time,end=end_time,freq='H').strftime('%Y%m%d%H')
            dfiles = []
            for j,dhour in enumerate(dhours):
                dfile = os.path.join(dlist,'wrfout_d01_'+dhour)
                if os.path.isfile(dfile):
    #                print(dfile)
                    dfiles.append(dfile)
            data_reinit_sn = xr.open_mfdataset(dfiles,combine='nested',concat_dim='Time')
            hours = np.array([float(hh) for hh in hours])
            times[count:count+len(dfiles)-1] = hours[1:len(dfiles)]
            VAR[count:count+len(dfiles)-1] = data_reinit_sn[var][1:len(dfiles)].values
        else :
            raise NameError
        count += 24
    root.close()
    
    #%%
    nc_file_path = dpath3 + 'WRF_reinit_9km_h_sfc_{}_{}.nc'.format(var_file,year)
    data = xr.open_dataset(nc_file_path)
    VAR_hourly = data[var].values
    dates = pd.date_range('{}-01-01'.format(year),'{}-12-31'.format(year),freq='24H').strftime('%Y%m%d%H')
    dates = [int(date) for date in dates]
    VAR_hourly = VAR_hourly.reshape(len(dates),24,360,530)
    
    nc_file_path = dpath3 + 'WRF_reinit_9km_d_sfc_{}_{}.nc'.format(var_file,year)
    root = nc.Dataset(nc_file_path,mode = 'w',format='NETCDF4')
    # 为nc文件创建维度
    lon = root.createDimension('lon',530)
    lat = root.createDimension('lat',360)
    time = root.createDimension('time',len(dates))
    
    # 为nc文件创建变量
    times = root.createVariable('times',int,('time'))
    longitude = root.createVariable('longitude',np.float32,('lat','lon'))
    latitude = root.createVariable('latitude',np.float32,('lat','lon'))
    
    VAR = root.createVariable(var,np.float32,('time','lat','lon'))
    # 对nc文件进行描述
    root.description = 'Nanjing University,Science of Atmosphere'
    root.history = 'Created ' + tm.ctime(tm.time())
    root.source = 'netCDF4 python '
    latitude.units = 'degree_nonth'
    longitude.units = 'degree_east'
    VAR.units = unit
    
    longitude[:,:] = lon_reinit_sn
    latitude[:,:] = lat_reinit_sn
    
    times[:len(dates)] = dates
    VAR[:len(dates)] = np.average(VAR_hourly,axis = 1)
    
    root.close()
    
    #%%
    nc_file_path = dpath3 + 'WRF_reinit_9km_d_sfc_{}_{}.nc'.format(var_file,year)
    data = xr.open_dataset(nc_file_path)
    VAR_daily = data[var].values
    dates = pd.date_range('{}-01-01'.format(year),'{}-12-31'.format(year),freq='24H').strftime('%Y%m%d')
    index_month = []
    for date in dates :
        if date[-2:] == '01':
            index_month.append(dates.get_loc(date))
    index_month.append(len(dates))
    
    
    nc_file_path = dpath3 + 'WRF_reinit_9km_m_sfc_{}_{}.nc'.format(var_file,year)
    root = nc.Dataset(nc_file_path,mode = 'w',format='NETCDF4')
    # 为nc文件创建维度
    lon = root.createDimension('lon',530)
    lat = root.createDimension('lat',360)
    time = root.createDimension('time',12)
    
    # 为nc文件创建变量
    times = root.createVariable('times',int,('time'))
    longitude = root.createVariable('longitude',np.float32,('lat','lon'))
    latitude = root.createVariable('latitude',np.float32,('lat','lon'))
    
    VAR = root.createVariable(var,np.float32,('time','lat','lon'))
    # 对nc文件进行描述
    root.description = 'Nanjing University,Science of Atmosphere'
    root.history = 'Created ' + tm.ctime(tm.time())
    root.source = 'netCDF4 python '
    latitude.units = 'degree_nonth'
    longitude.units = 'degree_east'
    VAR.units = unit
    
    longitude[:,:] = lon_reinit_sn
    latitude[:,:] = lat_reinit_sn
    for i in range(12): 
        times[i] = i+1
        VAR[i] = np.average(VAR_daily[index_month[i]:index_month[i+1]],axis = 0)
    root.close()
    
    #%%
    nc_file_path = dpath3 + 'WRF_reinit_9km_y_sfc_{}_{}.nc'.format(var_file,year)
    root = nc.Dataset(nc_file_path,mode = 'w',format='NETCDF4')
    # 为nc文件创建维度
    lon = root.createDimension('lon',530)
    lat = root.createDimension('lat',360)
    time = root.createDimension('time',1)
    
    # 为nc文件创建变量
    times = root.createVariable('times',int,('time'))
    longitude = root.createVariable('longitude',np.float32,('lat','lon'))
    latitude = root.createVariable('latitude',np.float32,('lat','lon'))
    
    VAR = root.createVariable(var,np.float32,('time','lat','lon'))
    # 对nc文件进行描述
    root.description = 'Nanjing University,Science of Atmosphere'
    root.history = 'Created ' + tm.ctime(tm.time())
    root.source = 'netCDF4 python '
    latitude.units = 'degree_nonth'
    longitude.units = 'degree_east'
    VAR.units = unit
    
    longitude[:,:] = lon_reinit_sn
    latitude[:,:] = lat_reinit_sn

    times[0] = year
    VAR[0] = np.average(VAR_daily,axis = 0)
    root.close()
