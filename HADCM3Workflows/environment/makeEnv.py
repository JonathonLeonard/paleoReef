import os
import pygmt
import argparse
import xarray as xr
import numpy as np
import xesmf as xe
import pandas as pd
import rioxarray as rio
from scipy.ndimage import gaussian_filter

def regriddingResolutionFine(res):
    ds_in = xr.DataArray(
        dims=["latitude", "longitude"],
        coords={"latitude": np.arange(-90,92.5,2.5), "longitude": np.arange(0,360,3.75)},
    )
    if res == 'res1':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.1, -90.0, 90.0, 0.1)
    elif res == 'res2':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
    elif res == 'res3':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.5, -90.0, 90.0, 0.5) 
    elif res == 'res4':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 1, -90.0, 90.0, 1) 
    elif res == 'res5':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 2, -90.0, 90.0, 2) 
    else:
        ds_out = xe.util.grid_2d(-180.0, 180.0, 3, -90.0, 90.0, 3)
    regrid = xe.Regridder(ds_in, ds_out, 'bilinear', periodic=True, weights='data/regrid/bilinear_gridder_'+res+'.nc')  
    return regrid

def regriddingResolutionCellFine(res):
    ds_in = xr.DataArray(
        dims=["latitude", "longitude"],
        coords={"latitude": np.arange(-88.75,91.25,2.5), "longitude": np.arange(1.875,361.875,3.75)},
    )
    if res == 'res1':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.1, -90.0, 90.0, 0.1)
    elif res == 'res2':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
    elif res == 'res3':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.5, -90.0, 90.0, 0.5)
    elif res == 'res4':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 1, -90.0, 90.0, 1)   
    elif res == 'res5':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 2, -90.0, 90.0, 2)   
    else:
        ds_out = xe.util.grid_2d(-180.0, 180.0, 3, -90.0, 90.0, 3)     
    regrid = xe.Regridder(ds_in, ds_out, 'bilinear', periodic=True, weights='data/regrid/cell_bilinear_gridder_'+res+'.nc') 
    return regrid

def regriddingInput(data, res):

    ds_in = xr.DataArray(
        dims=["latitude", "longitude"],
        coords={"latitude": np.arange(-90,90.05,0.05), "longitude": np.arange(-180,180.05,0.05)},
    )
    if res == 'res1':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.1, -90.0, 90.0, 0.1)
    elif res == 'res2':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
    elif res == 'res3':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.5, -90.0, 90.0, 0.5)
    elif res == 'res4':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 1, -90.0, 90.0, 1)
    elif res == 'res5':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 2, -90.0, 90.0, 2)
    else:
        ds_out = xe.util.grid_2d(-180.0, 180.0, 3, -90.0, 90.0, 3)   
    regrid = xe.Regridder(ds_in, ds_out, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_'+res+'.nc')
    rda = regrid(data)
      
    nda = xr.DataArray(
                rda.flood.values,
                dims=["latitude", "longitude"],
                coords={"latitude": rda.lat.values[:,0], "longitude": rda.lon.values[0,:]},
            )
    
    return nda.to_dataset(name='flood'), regrid

def regriddingResolutionCoarse(carbdata, maxdepth, res):
    ds_in = xr.DataArray(
        dims=["latitude", "longitude"],
        coords={"latitude": np.arange(-90,90.05,0.05), "longitude": np.arange(-180,180.05,0.05)},
    )
    if res == 'res1':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.1, -90.0, 90.0, 0.1)
    elif res == 'res2':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
    elif res == 'res3':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.5, -90.0, 90.0, 0.5)
    elif res == 'res4':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 1, -90.0, 90.0, 1)
    elif res == 'res5':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 2, -90.0, 90.0, 2)
    else:
        ds_out = xe.util.grid_2d(-180.0, 180.0, 3, -90.0, 90.0, 3)   
    regrid = xe.Regridder(ds_in, ds_out, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_'+res+'.nc')
    carb_da = regrid(carbdata)
    carb_da = carb_da.where(carb_da.flood<1.9) 
    carb_da = carb_da.where(carb_da.elevation>maxdepth)
    
    return carb_da

def regriddingResolutionCoarse2(carbdata, res):
    
    ds_in = xr.DataArray(
        dims=["latitude", "longitude"],
        coords={"latitude": np.arange(-90,90.05,0.05), "longitude": np.arange(-180,180.05,0.05)},
    )

    if res == 'res1':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.1, -90.0, 90.0, 0.1)
    elif res == 'res2':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
    elif res == 'res3':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.5, -90.0, 90.0, 0.5)
    elif res == 'res4':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 1, -90.0, 90.0, 1)
    elif res == 'res5':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 2, -90.0, 90.0, 2)
    else:
        ds_out = xe.util.grid_2d(-180.0, 180.0, 3, -90.0, 90.0, 3)  
    regrid = xe.Regridder(ds_in, ds_out, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_'+res+'.nc')
    
    return regrid(carbdata)

def getDataEnvi1(timeclim, month, res):
    
    ds_in = xr.DataArray(
        dims=["latitude", "longitude"],
        coords={"latitude": np.arange(-90,92.5,2.5), "longitude": np.arange(0,360,3.75)},
    )
    
    if res == 'res1':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.1, -90.0, 90.0, 0.1)
    elif res == 'res2':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
    elif res == 'res3':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.5, -90.0, 90.0, 0.5)
    elif res == 'res4':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 1, -90.0, 90.0, 1)
    elif res == 'res5':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 2, -90.0, 90.0, 2)
    else:
        ds_out = xe.util.grid_2d(-180.0, 180.0, 3, -90.0, 90.0, 3)  
    regrid = xe.Regridder(ds_in, ds_out, 'bilinear', periodic=True, weights='data/regrid/bilinear_gridder_'+res+'.nc')

    cds = xr.open_dataset('data/'+CO2curve+'/'+str(timeclim)+'Ma_'+months[month]+'.nc',decode_times=False).isel(t=0,depth_1=0,depth=0,unspecified=0)
    cds.rio.write_crs("epsg:4326", inplace=True)
    da = cds[envi_var1]
    da = da.rio.interpolate_na(method='cubic')
    da = da.rio.interpolate_na()
    rda = regrid(da)
    
    datalist = []
    for k in range(len(envi_var1)):
        if envi_var1[k] == 'salinity_mm_dpth':
            tmp = rda[envi_var1[k]].copy()
            tmp = tmp.where(tmp>-30/1000,-30/1000)
            datalist.append(tmp.values.copy())
        else:
            datalist.append(rda[envi_var1[k]].values.copy())
        
    return datalist

def getDataEnvi2(timeclim, month, res):
    
    cds = xr.open_dataset('data/'+CO2curve+'/'+str(timeclim)+'Ma_'+months[month]+'.nc',decode_times=False).isel(t=0,depth_1=0,depth=0,unspecified=0)
    
    da = cds[envi_var2] 
    cur = np.sqrt(da[envi_var2[0]].values**2+da[envi_var2[1]].values**2)
    dcur = xr.DataArray(
                cur,
                dims=["latitude", "longitude"],
                coords={"latitude": cds.latitude_1.values, "longitude": cds.longitude_1.values},
            )
    dcur.rio.write_crs("epsg:4326", inplace=True)
    dcur.rio.write_nodata(np.nan, inplace=True)
    dcur_interp = dcur.rio.interpolate_na(method='cubic')
    dcur_interp = dcur_interp.rio.interpolate_na()
    dcur_interp = dcur_interp.where(dcur_interp>0,0)
    ds_in = xr.DataArray(
            dims=["latitude", "longitude"],
            coords={"latitude": np.arange(-88.75,91.25,2.5), "longitude": np.arange(1.875,361.875,3.75)},
        )
    if res == 'res1':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.1, -90.0, 90.0, 0.1)
    elif res == 'res2':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
    elif res == 'res3':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 0.5, -90.0, 90.0, 0.5)
    elif res == 'res4':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 1, -90.0, 90.0, 1)
    elif res == 'res5':
        ds_out = xe.util.grid_2d(-180.0, 180.0, 2, -90.0, 90.0, 2)
    else:
        ds_out = xe.util.grid_2d(-180.0, 180.0, 3, -90.0, 90.0, 3)  
    regrid = xe.Regridder(ds_in, ds_out, 'bilinear', periodic=True, weights='data/regrid/cell_bilinear_gridder_'+res+'.nc')

    return regrid(dcur_interp).values.copy()

def saveGeoTiff(cgeo,var,reso):
    gtmp = cgeo[var].rio.write_nodata(np.nan, inplace=True).to_dataset()
    gtmp = xr.DataArray(
                    gtmp[var].values,
                    dims=["latitude", "longitude"],
                    coords={"latitude": gtmp.lat.values[:,0], "longitude": gtmp.lon.values[0,:]},
                )
    gtmp = gtmp.to_dataset(name=var)
    gtmp.rio.write_crs('EPSG:4326', inplace=True)
    gtmp[var].rio.to_raster('env_var/bio_'+reso+'_'+str(timeclim)+'Ma_'+CO2curve+'/'+var+'.tif')

    return

################################################################
################################################################ 

# Parsing command line arguments
parser = argparse.ArgumentParser(
    description="This is a simple entry to adjust elevation.", add_help=True
)
parser.add_argument("-t", "--time", help="Time interval", required=True)
parser.add_argument("-c", "--co", help="CO2 curve", required=True)
# parser.add_argument("-r", "--res", help="Resolution grid", default='3')

args = parser.parse_args()

timeclim = int(args.time)
coarg = args.co
res = 'res3' #+args.res

CO2curve = 'smooth'
if coarg == '1':
    CO2curve = 'foster'

# Carbonate factory
dsCarbh = xr.open_dataset('data/carbonate_factory.nc')
# Shelf from Scotese
shelf = xr.open_dataset('data/shelfs/'+str(timeclim)+'Ma.nc')

# if res == 'res1':
#     print('Environmental variables grid resolution 0.1 degrees')
# elif res == 'res2':
#     print('Grid resolution 0.25 degrees')
# elif res == 'res3':
#     print('Grid resolution 0.5 degrees')
# elif res == 'res4':
#     print('Grid resolution 1.0 degrees')
# elif res == 'res5':
#     print('Grid resolution 2.0 degrees')
# elif res == 'res6':
#     print('Grid resolution 3.0 degrees')
# else:
#     print('Choice for resolution should be either [1,2,3,4,5]')
#     exit()

maxdepth = -50000.
varnames = ['photozoan'] #['biochemical', 'photozoan', 'photoC', 'heterozoan']
dsCarbh['flood'] = (('lat', 'lon'), shelf.flood.values)
shelf.close()

months = ['jan','feb','mar','apr','may','jun','jul','aug','sep','oct','nov','dec']
envi_var1 = ['temp_mm_dpth','salinity_mm_dpth','streamFn_mm_uo'] # 'mixLyrDpth_mm_uo','iceconc_mm_uo'
envi_var2 = ['ucurrTot_mm_dpth','vcurrTot_mm_dpth']

envi_var = ['temp','salt','fstream','curr'] # 'mixlay','ice'
# envi_attr = ['temp_mm_dpth','salinity_mm_dpth','streamFn_mm_uo','ucurrTot_mm_dpth'] # 'mixLyrDpth_mm_uo','iceconc_mm_uo'

regridFiner = regriddingResolutionFine(res)
regridFinerCell = regriddingResolutionCellFine(res)

dsCarb = dsCarbh.rename({'lon': 'longitude','lat': 'latitude'})
CarbFac = regriddingResolutionCoarse(dsCarb, maxdepth, res)

## Store each factory position
pathf = 'env_var/factoryReef'
isExist = os.path.exists(pathf)
if not isExist:
    os.makedirs(pathf)

# if timeclim == 0:
#     for k in range(len(varnames)):
#         tmp = CarbFac[varnames[k]].copy()
#         lon = tmp.lon.values[0,:] 
#         lat = tmp.lat.values[:,0]
#         ids = np.stack(np.argwhere(~np.isnan(tmp.values)),axis=0)
#         coords = np.zeros(ids.shape)
#         coords[:,0] =  np.round(lon[ids[:,1]],2)
#         coords[:,1] =  np.round(lat[ids[:,0]],2)
#         data = {'lon': coords[:,0],
#                 'lat': coords[:,1]}
#         df = pd.DataFrame(data)
        # df.to_csv(pathf+'/'+varnames[k]+'_'+res+'.csv',index=False)

env1_monthly = [] 
env2_monthly = [] 
for k in range(12):
    env1_monthly.append(getDataEnvi1(timeclim, k, res)) 
    env2_monthly.append(getDataEnvi2(timeclim, k, res)) 

minvar = []
maxvar = []
meanvar = []
for v in range(len(envi_var1)):
    tmp = env1_monthly[0][v].copy()
    tmp_min = env1_monthly[0][v].copy()
    tmp_max = env1_monthly[0][v].copy()
    for k in range(1,12):
        tmp_max = np.maximum(tmp_max,env1_monthly[k][v])
        tmp_min = np.minimum(tmp_min,env1_monthly[k][v])
        tmp += env1_monthly[k][v]
    maxvar.append(tmp_max)  
    minvar.append(tmp_min)  
    meanvar.append(tmp/12.)  
    
tmp_max = env2_monthly[0]
tmp_min = env2_monthly[0]
tmp = env2_monthly[0]
for k in range(1,12):
    tmp_max = np.maximum(tmp_max,env2_monthly[k])
    tmp_min = np.minimum(tmp_min,env2_monthly[k])
    tmp += env2_monthly[k]
maxvar.append(tmp_max)  
minvar.append(tmp_min)  
meanvar.append(tmp/12.) 

cds = xr.open_dataset('data/'+CO2curve+'/'+str(timeclim)+'Ma_'+months[0]+'.nc',decode_times=False).isel(t=0,depth_1=0,depth=0,unspecified=0)
abiotic, regrid = regriddingInput(xr.open_dataset('data/shelfs/'+str(timeclim)+'Ma.nc'),res)

if timeclim == 0:
    tmpda = xr.open_dataset('data/carbonate_factory.nc')
    tmpda = tmpda.rename({'lon': 'longitude','lat': 'latitude'})
    bathy = regrid(tmpda.elevation)
    bathy = bathy.where(bathy>=0,0)
    abiotic['bathy'] = (('latitude', 'longitude'), bathy.values)
else:
    tmpda = xr.open_dataset('data/physio/physio'+str(timeclim)+'Ma.nc') 
    bathy = regrid(tmpda.zpaleo)
    bathy = bathy.where(bathy>=0,0)
    abiotic['bathy'] = (('latitude', 'longitude'), bathy.values)
abiotic['bathy'].attrs['name'] = 'z'

for k in range(len(envi_var)):
    fac1 = 1
    fac2 = 0
    if envi_var[k] == 'salt':
        fac1 = 1000
        fac2 = 35
    if envi_var[k] == 'fstream':
        fac1 = 1e-12
    if envi_var[k] == 'ice':
        fac1 = 100
    lbd = 20
    abiotic[envi_var[k]+'Min'] = (('latitude', 'longitude'), gaussian_filter(minvar[k]*fac1+fac2,lbd)) #minvar[k]*fac1+fac2
    abiotic[envi_var[k]+'Min'].attrs['name'] = envi_var[k]+'min'
    abiotic[envi_var[k]+'Max'] = (('latitude', 'longitude'), gaussian_filter(maxvar[k]*fac1+fac2,lbd)) #maxvar[k]*fac1+fac2)
    abiotic[envi_var[k]+'Max'].attrs['name'] = envi_var[k]+'max'
    abiotic[envi_var[k]] = (('latitude', 'longitude'), gaussian_filter(meanvar[k]*fac1+fac2,lbd)) #meanvar[k]*fac1+fac2)
    abiotic[envi_var[k]].attrs['name'] = envi_var[k]

gridLat = abiotic.latitude.values.copy()
gridLon = abiotic.longitude.values.copy()
seddf = pd.read_csv('data/wsflux/sed'+str(timeclim)+'.csv')
sed = seddf['val'].values*1.e9 # m3 /yr
logSed = np.log10(sed)
sedLon = seddf['lon'].values
sedLat = seddf['lat'].values

sedflx = np.zeros(bathy.values.shape)
for k in range(len(seddf)):
    lonlat = abiotic.sel(longitude=sedLon[k],latitude=sedLat[k],method='nearest')
    idlat = np.argwhere(gridLat==lonlat.latitude.values)
    idlon = np.argwhere(gridLon==lonlat.longitude.values)
    sedflx[idlat,idlon] = sed[k]

lbd1 = 80
lbd2 = 10
sedflxblur = gaussian_filter(sedflx, sigma=lbd1) 
abiotic['sed'] = (('latitude', 'longitude'), sedflxblur)
val = abiotic.sed.where(abiotic.bathy>maxdepth,0)
newsed = val.values.copy()
newsed = gaussian_filter(newsed, sigma=lbd2) 
abiotic['sed'] = (('latitude', 'longitude'), newsed)
val = abiotic.sed.where(abiotic.flood<1.2)
abiotic['sed'] = val
abiotic['sed'].attrs['name'] = 'river'

abiotic = abiotic.where(abiotic.flood<1.9)
try:
    os.remove('env_var/enviVar'+str(timeclim)+'Ma_'+res+'_'+CO2curve+'.nc')
except FileNotFoundError:
    print("Save netCDF file.")

abiotic.to_netcdf('env_var/enviVar'+str(timeclim)+'Ma_'+res+'_'+CO2curve+'.nc')

geoda = abiotic.copy()
geoda = geoda.where(geoda.bathy>maxdepth)

# Make coarser (0.25 degree) environmental outputs
ds_in = xr.DataArray(
        dims=["latitude", "longitude"],
        coords={"latitude": np.arange(-89.95,90.05,0.1), "longitude": np.arange(-179.95,180.05,0.1)},
    )
ds_out2 = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
ds_out3 = xe.util.grid_2d(-180.0, 180.0, 0.5, -90.0, 90.0, 0.5)
ds_out4 = xe.util.grid_2d(-180.0, 180.0, 1, -90.0, 90.0, 1)
ds_out5 = xe.util.grid_2d(-180.0, 180.0, 2, -90.0, 90.0, 2)
ds_out6 = xe.util.grid_2d(-180.0, 180.0, 3, -90.0, 90.0, 3)
regrid2 = xe.Regridder(ds_in, ds_out2, 'bilinear', periodic=True, weights='data/regrid/bilinear_gridder_res2.nc')
regrid3 = xe.Regridder(ds_in, ds_out3, 'bilinear', periodic=True, weights='data/regrid/bilinear_gridder_res3.nc')
regrid4 = xe.Regridder(ds_in, ds_out4, 'bilinear', periodic=True, weights='data/regrid/bilinear_gridder_res4.nc')
regrid5 = xe.Regridder(ds_in, ds_out5, 'bilinear', periodic=True, weights='data/regrid/bilinear_gridder_res5.nc')
regrid6 = xe.Regridder(ds_in, ds_out6, 'bilinear', periodic=True, weights='data/regrid/bilinear_gridder_res6.nc')

regridn2 = xe.Regridder(ds_in, ds_out2, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_res2.nc')
regridn3 = xe.Regridder(ds_in, ds_out3, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_res3.nc')
regridn4 = xe.Regridder(ds_in, ds_out4, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_res4.nc')
regridn5 = xe.Regridder(ds_in, ds_out5, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_res5.nc')
regridn6 = xe.Regridder(ds_in, ds_out6, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_res6.nc')


# Remove all data not on shelf
# geoda = geoda.where(geoda.shelf==1)
coarse_geoda2 = regrid2(geoda)
coarse_geoda3 = regrid3(geoda)
coarse_geoda4 = regrid4(geoda)
coarse_geoda5 = regrid5(geoda)
coarse_geoda6 = regrid6(geoda)

path = 'env_var/bio_res1_'+str(timeclim)+'Ma_'+CO2curve
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
path = 'env_var/bio_res2_'+str(timeclim)+'Ma_'+CO2curve
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
path = 'env_var/bio_res3_'+str(timeclim)+'Ma_'+CO2curve
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
path = 'env_var/bio_res4_'+str(timeclim)+'Ma_'+CO2curve
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
path = 'env_var/bio_res5_'+str(timeclim)+'Ma_'+CO2curve
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)
path = 'env_var/bio_res6_'+str(timeclim)+'Ma_'+CO2curve
isExist = os.path.exists(path)
if not isExist:
    os.makedirs(path)

# Bathymetry
geotmp = geoda['bathy'].rio.write_nodata(np.nan, inplace=True).to_dataset()
geotmp.rio.write_crs('EPSG:4326', inplace=True)
geotmp['bathy'].rio.to_raster('env_var/bio_res1_'+str(timeclim)+'Ma_'+CO2curve+'/bathy.tif')
saveGeoTiff(coarse_geoda2,'bathy','res2')
saveGeoTiff(coarse_geoda3,'bathy','res3')
saveGeoTiff(coarse_geoda4,'bathy','res4')
saveGeoTiff(coarse_geoda5,'bathy','res5')
saveGeoTiff(coarse_geoda6,'bathy','res6')

# Terrigenous flux
geotmp = geoda['sed'].rio.write_nodata(np.nan, inplace=True).to_dataset()
geotmp.rio.write_crs('EPSG:4326', inplace=True)
geotmp['sed'].rio.to_raster('env_var/bio_res1_'+str(timeclim)+'Ma_'+CO2curve+'/sed.tif')
saveGeoTiff(coarse_geoda2,'sed','res2')
saveGeoTiff(coarse_geoda3,'sed','res3')
saveGeoTiff(coarse_geoda4,'sed','res4')
saveGeoTiff(coarse_geoda5,'sed','res5')
saveGeoTiff(coarse_geoda6,'sed','res6')

# Other variables
envi_var = ['temp','salt','fstream','curr']
for k in range(len(envi_var)):
    if envi_var[k] == 'temp' or envi_var[k] == 'salt':
        geotmp = geoda[envi_var[k]+'Min'].rio.write_nodata(np.nan, inplace=True).to_dataset()
        geotmp.rio.write_crs('EPSG:4326', inplace=True)
        geotmp[envi_var[k]+'Min'].rio.to_raster('env_var/bio_res1_'+str(timeclim)+'Ma_'+CO2curve+'/'+envi_var[k]+"Min.tif") 
        saveGeoTiff(coarse_geoda2,envi_var[k]+'Min','res2')
        saveGeoTiff(coarse_geoda3,envi_var[k]+'Min','res3')
        saveGeoTiff(coarse_geoda4,envi_var[k]+'Min','res4')
        saveGeoTiff(coarse_geoda5,envi_var[k]+'Min','res5')
        saveGeoTiff(coarse_geoda6,envi_var[k]+'Min','res6')

        geotmp = geoda[envi_var[k]+'Max'].rio.write_nodata(np.nan, inplace=True).to_dataset()
        geotmp.rio.write_crs('EPSG:4326', inplace=True)
        geotmp[envi_var[k]+'Max'].rio.to_raster('env_var/bio_res1_'+str(timeclim)+'Ma_'+CO2curve+'/'+envi_var[k]+"Max.tif")
        saveGeoTiff(coarse_geoda2,envi_var[k]+'Max','res2')
        saveGeoTiff(coarse_geoda3,envi_var[k]+'Max','res3')
        saveGeoTiff(coarse_geoda4,envi_var[k]+'Max','res4')
        saveGeoTiff(coarse_geoda5,envi_var[k]+'Max','res5')
        saveGeoTiff(coarse_geoda6,envi_var[k]+'Max','res6')
    else:
        geotmp = geoda[envi_var[k]].rio.write_nodata(np.nan, inplace=True).to_dataset()
        geotmp.rio.write_crs('EPSG:4326', inplace=True)
        geotmp[envi_var[k]].rio.to_raster('env_var/bio_res1_'+str(timeclim)+'Ma_'+CO2curve+'/'+envi_var[k]+".tif")
        saveGeoTiff(coarse_geoda2,envi_var[k],'res2')
        saveGeoTiff(coarse_geoda3,envi_var[k],'res3')
        saveGeoTiff(coarse_geoda4,envi_var[k],'res4')
        saveGeoTiff(coarse_geoda5,envi_var[k],'res5')
        saveGeoTiff(coarse_geoda6,envi_var[k],'res6')

if timeclim == 0:
    dsCarbh = xr.open_dataset('data/carbonate_factory.nc')
    shelf = xr.open_dataset('data/shelfs/'+str(timeclim)+'Ma.nc')
    dsCarbh['flood'] = (('lat', 'lon'), shelf.flood.values)
    shelf.close()
    dsCarb = dsCarbh.rename({'lon': 'longitude','lat': 'latitude'})
    CarbFac = regriddingResolutionCoarse2(dsCarb, res)
    abiotic['elevation'] = (('latitude', 'longitude'), CarbFac.elevation.values)
    abiotic['photozoan'] = (('latitude', 'longitude'), CarbFac.photozoan.values)
    abiotic['photoC'] = (('latitude', 'longitude'), CarbFac.photoC.values)
    abiotic['heterozoan'] = (('latitude', 'longitude'), CarbFac.heterozoan.values)
    abiotic['biochemical'] = (('latitude', 'longitude'), CarbFac.biochemical.values)
    abiotic['shelf'] = (('latitude', 'longitude'), CarbFac.flood.values)

    val = abiotic[['photozoan','photoC','heterozoan','biochemical','elevation']].fillna(0)
    val = val.where(val <= 0, 1).where(val.elevation<0)
    tmp = val.where(val == 0, 0).where(val.elevation>-2500).photozoan
    tmp1 = val.photozoan.where((val.elevation>-1000)&(val.photozoan>0)&(val.elevation != np.nan)).fillna(0)
    (tmp+tmp1).to_dataframe().dropna().to_csv(pathf+'/res1_photozoan.csv')
    # val1.heterozoan.to_dataframe().dropna().to_csv(pathf+'/res1_heterozoan.csv')
    # val1.photoC.to_dataframe().dropna().to_csv(pathf+'/res1_photoC.csv')
    # val1.biochemical.to_dataframe().dropna().to_csv(pathf+'/res1_biochemical.csv')

    val = regridn2(abiotic[['photozoan','photoC','heterozoan','biochemical','elevation']].fillna(0))
    val = val.where(val <= 0, 1).where(val.elevation<0)
    tmp = val.where(val == 0, 0).where(val.elevation>-2500).photozoan
    tmp1 = val.photozoan.where((val.elevation>-1000)&(val.photozoan>0)&(val.elevation != np.nan)).fillna(0)
    (tmp+tmp1).to_dataframe().dropna().to_csv(pathf+'/res2_photozoan.csv')
    # val2.heterozoan.to_dataframe().dropna().to_csv(pathf+'/res2_heterozoan.csv')
    # val2.photoC.to_dataframe().dropna().to_csv(pathf+'/res2_photoC.csv')
    # val2.biochemical.to_dataframe().dropna().to_csv(pathf+'/res2_biochemical.csv')

    val = regridn3(abiotic[['photozoan','photoC','heterozoan','biochemical','elevation']].fillna(0))
    val = val.where(val <= 0, 1).where(val.elevation<0)
    tmp = val.where(val == 0, 0).where(val.elevation>-2500).photozoan
    tmp1 = val.photozoan.where((val.elevation>-1000)&(val.photozoan>0)&(val.elevation != np.nan)).fillna(0)
    (tmp+tmp1).to_dataframe().dropna().to_csv(pathf+'/res3_photozoan.csv')
    # val3.heterozoan.to_dataframe().dropna().to_csv(pathf+'/res3_heterozoan.csv')
    # val3.photoC.to_dataframe().dropna().to_csv(pathf+'/res3_photoC.csv')
    # val3.biochemical.to_dataframe().dropna().to_csv(pathf+'/res3_biochemical.csv')

    val = regridn4(abiotic[['photozoan','photoC','heterozoan','biochemical','elevation']].fillna(0))
    val = val.where(val <= 0, 1).where(val.elevation<0)
    tmp = val.where(val == 0, 0).where(val.elevation>-2500).photozoan
    tmp1 = val.photozoan.where((val.elevation>-1000)&(val.photozoan>0)&(val.elevation != np.nan)).fillna(0)
    (tmp+tmp1).to_dataframe().dropna().to_csv(pathf+'/res4_photozoan.csv')
    # val4.heterozoan.to_dataframe().dropna().to_csv(pathf+'/res4_heterozoan.csv')
    # val4.photoC.to_dataframe().dropna().to_csv(pathf+'/res4_photoC.csv')
    # val4.biochemical.to_dataframe().dropna().to_csv(pathf+'/res4_biochemical.csv')

    val = regridn5(abiotic[['photozoan','photoC','heterozoan','biochemical','elevation']].fillna(0))
    val = val.where(val <= 0, 1).where(val.elevation<0)
    tmp = val.where(val == 0, 0).where(val.elevation>-2500).photozoan
    tmp1 = val.photozoan.where((val.elevation>-1000)&(val.photozoan>0)&(val.elevation != np.nan)).fillna(0)
    (tmp+tmp1).to_dataframe().dropna().to_csv(pathf+'/res5_photozoan.csv')
    # val5.heterozoan.to_dataframe().dropna().to_csv(pathf+'/res5_heterozoan.csv')
    # val5.photoC.to_dataframe().dropna().to_csv(pathf+'/res5_photoC.csv')
    # val5.biochemical.to_dataframe().dropna().to_csv(pathf+'/res5_biochemical.csv')

    val = regridn6(abiotic[['photozoan','photoC','heterozoan','biochemical','elevation']].fillna(0))
    val = val.where(val <= 0, 1).where(val.elevation<0)
    tmp = val.where(val == 0, 0).where(val.elevation>-2500).photozoan
    tmp1 = val.photozoan.where((val.elevation>-1000)&(val.photozoan>0)&(val.elevation != np.nan)).fillna(0)
    (tmp+tmp1).to_dataframe().dropna().to_csv(pathf+'/res6_photozoan.csv')
    # val6.heterozoan.to_dataframe().dropna().to_csv(pathf+'/res6_heterozoan.csv')
    # val6.photoC.to_dataframe().dropna().to_csv(pathf+'/res6_photoC.csv')
    # val6.biochemical.to_dataframe().dropna().to_csv(pathf+'/res6_biochemical.csv')
else:
    abiotic['elevation'] = abiotic.bathy.copy()
    shelf = xr.open_dataset('data/shelfs/'+str(timeclim)+'Ma.nc')
    shelf2 = regriddingResolutionCoarse2(shelf, res)
    shelf.close()
    abiotic['shelf'] = (('latitude', 'longitude'), shelf2.flood.values)

# Plots
continent = abiotic.flood.fillna(10)
continent = continent.where(continent>4)
region = [-180,180,-80,80]#'d'

path2 = 'figs/env_'+res+'_'+str(timeclim)+'Ma_'+CO2curve
isExist = os.path.exists(path2)
if not isExist:
    os.makedirs(path2)

pygmt.config(
    MAP_FRAME_WIDTH="0.1p",
    MAP_FRAME_TYPE="fancy",
)
font = "4p,Helvetica-Bold"

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="panoply", series=[-5, 40])
    fig.grdimage(abiotic.temp, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a10", "x", r"y+l\260C"],
        scale=1,
        transparency=20,
    )    
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/temp_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="panoply", series=[-5, 40])
    fig.grdimage(abiotic.tempMax, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a10", "x", r"y+l\260C"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/tempMax_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)    
    pygmt.makecpt(cmap="panoply", series=[-5, 40])
    fig.grdimage(abiotic.tempMin, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a10", "x", r"y+l\260C"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/tempMin_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')    
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="panoply", series=[28, 40])
    fig.grdimage(abiotic.salt, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a2", "x", r"y+lPSU"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/salt_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')    
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="panoply", series=[28, 40])
    fig.grdimage(abiotic.saltMax, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a2", "x", r"y+lPSU"],
        scale=1,
        transparency=20,
    )    
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/saltMax_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')    
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="panoply", series=[28, 40])
    fig.grdimage(abiotic.saltMin, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a2", "x", r"y+lPSU"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/saltMin_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

# fig = pygmt.Figure()
# with pygmt.config(FONT='4p,Helvetica,black'):
#     pygmt.makecpt(cmap="gray", series=[-15000, 4000])
#     fig.basemap(region=region, projection="Q12c", frame='af')     
#     fig.grdimage(abiotic.elevation, shading=True, frame=False) 
#     pygmt.makecpt(cmap="vik", series=[0, 200, 5])
#     fig.grdimage(abiotic.mixlayMax, shading=False, nan_transparent=True, frame=False)
#     fig.colorbar(
#         position="g-170/-58+w2c/0.15c+h",
#         box="+gwhite+p0.01p",
#         frame=["a50", "x", r"y+lm"],
#         scale=1,
#         transparency=20,
#     )
#     pygmt.makecpt(cmap="black", series=[8, 10])
#     fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
#     fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
# fig.savefig(path2+'/mixlayMax_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

# fig = pygmt.Figure()
# with pygmt.config(FONT='4p,Helvetica,black'):
#     pygmt.makecpt(cmap="gray", series=[-15000, 4000])
#     fig.basemap(region=region, projection="Q12c", frame='af')     
#     fig.grdimage(abiotic.elevation, shading=True, frame=False) 
#     pygmt.makecpt(cmap="vik", series=[0, 200, 5])
#     fig.grdimage(abiotic.mixlayMin, shading=False, nan_transparent=True, frame=False)
#     fig.colorbar(
#         position="g-170/-58+w2c/0.15c+h",
#         box="+gwhite+p0.01p",
#         frame=["a50", "x", r"y+lm"],
#         scale=1,
#         transparency=20,
#     )
#     pygmt.makecpt(cmap="black", series=[8, 10])
#     fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
#     fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
# fig.savefig(path2+'/mixlayMin_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

# fig = pygmt.Figure()
# with pygmt.config(FONT='4p,Helvetica,black'):
#     pygmt.makecpt(cmap="gray", series=[-15000, 4000])
#     fig.basemap(region=region, projection="Q12c", frame='af')     
#     fig.grdimage(abiotic.elevation, shading=True, frame=False) 
#     pygmt.makecpt(cmap="vik", series=[0, 200, 5])
#     fig.grdimage(abiotic.mixlayMean, shading=False, nan_transparent=True, frame=False)
#     fig.colorbar(
#         position="g-170/-58+w2c/0.15c+h",
#         box="+gwhite+p0.01p",
#         frame=["a50", "x", r"y+lm"],
#         scale=1,
#         transparency=20,
#     )
#     pygmt.makecpt(cmap="black", series=[8, 10])
#     fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
#     fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
# fig.savefig(path2+'/mixlayMean_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')    
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="imola", series=[0, 30, 5])
    fig.grdimage(abiotic.curr, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a5", "x", r"y+lcm/s"],
        scale=1,
        transparency=20,
    )    
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/curr_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')    
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="imola", series=[0, 30, 5])
    fig.grdimage(abiotic.currMax, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a5", "x", r"y+lcm/s"],
        scale=1,
        transparency=20,
    )    
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/currMax_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')    
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="imola", series=[0, 30, 5])
    fig.grdimage(abiotic.currMin, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a5", "x", r"y+lcm/s"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/currMin_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

# fig = pygmt.Figure()
# with pygmt.config(FONT='4p,Helvetica,black'):
#     pygmt.makecpt(cmap="gray", series=[-15000, 4000])
#     fig.basemap(region=region, projection="Q12c", frame='af')     
#     fig.grdimage(abiotic.elevation, shading=True, frame=False) 
#     pygmt.makecpt(cmap="oslo", series=[0, 100], reverse=True)
#     val = abiotic.iceMax.copy()
#     val = val.where(val>0.1)
#     fig.grdimage(val, shading=False, nan_transparent=True, frame=False)
#     fig.colorbar(
#         position="g-170/-58+w2c/0.15c+h",
#         box="+gwhite+p0.01p",
#         frame=["a25", "x", r"y+l%"],
#         scale=1,
#         transparency=20,
#     )
#     pygmt.makecpt(cmap="black", series=[8, 10])
#     fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
#     fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
# fig.savefig(path2+'/iceMax_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

# fig = pygmt.Figure()
# with pygmt.config(FONT='4p,Helvetica,black'):
#     pygmt.makecpt(cmap="gray", series=[-15000, 4000])
#     fig.basemap(region=region, projection="Q12c", frame='af')     
#     fig.grdimage(abiotic.elevation, shading=True, frame=False) 
#     pygmt.makecpt(cmap="oslo", series=[0, 100], reverse=True)
#     val = abiotic.iceMin.copy()
#     val = val.where(val>0.1)
#     fig.grdimage(val, shading=False, nan_transparent=True, frame=False)
#     fig.colorbar(
#         position="g-170/-58+w2c/0.15c+h",
#         box="+gwhite+p0.01p",
#         frame=["a25", "x", r"y+l%"],
#         scale=1,
#         transparency=20,
#     )
#     pygmt.makecpt(cmap="black", series=[8, 10])
#     fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
#     fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
# fig.savefig(path2+'/iceMin_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

# fig = pygmt.Figure()
# with pygmt.config(FONT='4p,Helvetica,black'):
#     pygmt.makecpt(cmap="gray", series=[-15000, 4000])
#     fig.basemap(region=region, projection="Q12c", frame='af')     
#     fig.grdimage(abiotic.elevation, shading=True, frame=False) 
#     pygmt.makecpt(cmap="oslo", series=[0, 100], reverse=True)
#     val = abiotic.iceMean.copy()
#     val = val.where(val>0.1)
#     fig.grdimage(val, shading=False, nan_transparent=True, frame=False)
#     fig.colorbar(
#         position="g-170/-58+w2c/0.15c+h",
#         box="+gwhite+p0.01p",
#         frame=["a25", "x", r"y+l%"],
#         scale=1,
#         transparency=20,
#     )
#     pygmt.makecpt(cmap="black", series=[8, 10])
#     fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
#     fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
# fig.savefig(path2+'/iceMean_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af')     
    fig.grdimage(abiotic.elevation, shading=True, frame=False) 
    pygmt.makecpt(cmap="bam", series=[-200, 200, 25])
    fig.grdimage(abiotic.fstream, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a100", "x", r"y+lSv"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/fstream_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af')     
    fig.grdimage(abiotic.elevation, shading=True, frame=False) 
    pygmt.makecpt(cmap="bam", series=[-200, 200, 25])
    fig.grdimage(abiotic.fstreamMax, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a100", "x", r"y+lSv"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/fstreamMax_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af')     
    fig.grdimage(abiotic.elevation, shading=True, frame=False) 
    pygmt.makecpt(cmap="bam", series=[-200, 200, 25])
    fig.grdimage(abiotic.fstreamMin, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a100", "x", r"y+lSv"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[8, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/fstreamMin_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

seddf = pd.read_csv('data/wsflux/sed'+str(timeclim)+'.csv')
sed = seddf['val'].values*1.e9 # m3 /yr
logSed = np.log10(sed)
sorted_index_array1 = np.argsort(logSed)
sortedSed = logSed[sorted_index_array1]
sortedLon1 = seddf['lon'].values[sorted_index_array1]
sortedLat1 = seddf['lat'].values[sorted_index_array1]
nlargest = 500
rLon1 = sortedLon1[-nlargest : ]
rLat1 = sortedLat1[-nlargest : ]
rSed = sortedSed[-nlargest : ]
sedriver = abiotic.sed.where((abiotic.sed>10) & (abiotic.bathy>-1000))

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="gray", series=[-15000, 4000])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    pygmt.makecpt(cmap="davos", series=[0, 5.e3, 500], reverse=True)
    fig.grdimage(sedriver, shading=False, nan_transparent=True, frame=False)
    fig.colorbar(
        position="g-170/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a1000", "x", r"y+lkm3/yr"],
        scale=1,
        transparency=20,
    )
    pygmt.makecpt(cmap="black", series=[0, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
with pygmt.config(FONT='6p,Helvetica,black'):
    fig.plot(
        x=rLon1,
        y=rLat1,
        style="cc",
        pen="black",
        size=0.00075 * 2 ** rSed,
        fill="white",
    )
fig.savefig(path2+'/rivers_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)

fig = pygmt.Figure()
with pygmt.config(FONT='4p,Helvetica,black'):
    pygmt.makecpt(cmap="ibcso", series=[-10000, 0])
    fig.basemap(region=region, projection="Q12c", frame='af') #'afg')
    fig.grdimage(abiotic.elevation, shading=True, frame=False) #,transparency=50)
    fig.colorbar(
        position="g-165/-58+w2c/0.15c+h",
        box="+gwhite+p0.01p",
        frame=["a5000", "x", r"y+lm"],
        scale=1,
        transparency=20,
    )    
    pygmt.makecpt(cmap="black", series=[0, 10])
    fig.grdimage(continent, shading=True, nan_transparent=True, frame=False)
    fig.grdcontour(interval=1,grid=abiotic.shelf,limit=[1.9, 2],pen="0.5p,white")
fig.savefig(path2+'/bathy_'+str(timeclim)+'Ma_'+CO2curve+'.pdf',dpi=300)