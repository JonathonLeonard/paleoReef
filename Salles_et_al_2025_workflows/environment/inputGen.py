import os
import argparse
import numpy as np
import xesmf as xe
import xarray as xr
import pandas as pd
import rioxarray as rio


# Parsing command line arguments
parser = argparse.ArgumentParser(
    description="This is a simple entry to build env. data for biomod.", add_help=True
)
# parser.add_argument("-t", "--time", help="Time interval", required=True)
parser.add_argument("-c", "--co", help="CO2 curve", required=True)

args = parser.parse_args()

# timeclim = int(args.time)
coarg = args.co

CO2curve = 'smooth'
if coarg == '1':
    CO2curve = 'foster'

times = pd.read_csv('data/bristol_sim_list.csv')['time'].values.tolist()

for t in range(len(times)):
    timeclim = times[t]
    abiotic = xr.open_dataset('env_var/enviVar'+str(timeclim)+'Ma_res1_'+CO2curve+'.nc')
    flood = abiotic.flood.values.copy()
    bathy = abiotic.bathy.values.copy()

    ids = np.where((flood==1) & (bathy>=0))
    bathy[ids] = -0.01
    abiotic['bathy'] = (('latitude', 'longitude'), bathy)

    # Regridding
    ds_out = xe.util.grid_2d(-180.0, 180.0, 0.25, -90.0, 90.0, 0.25)
    regrid = xe.Regridder(abiotic, ds_out, 'nearest_s2d', periodic=True, weights='data/regrid/nearest_s2d_gridder_env_res2.nc')

    abiotic2 = regrid(abiotic)
    shelfdata = abiotic2.where(abiotic2.flood==1)

    pathf = 'env_var/bio_res2_'+str(timeclim)+'Ma_'+CO2curve
    isExist = os.path.exists(pathf)
    if not isExist:
        os.makedirs(pathf)
        
        
    gtmp = shelfdata['sed'].rio.write_nodata(np.nan, inplace=True).to_dataset()
    gtmp = xr.DataArray(
                    gtmp['sed'].values,
                    dims=["latitude", "longitude"],
                    coords={"latitude": gtmp.lat.values[:,0], "longitude": gtmp.lon.values[0,:]},
                )
    gtmp = gtmp.to_dataset(name='sed')
    gtmp.rio.write_crs('EPSG:4326', inplace=True)
    gtmp['sed'].rio.to_raster('env_var/bio_res2_'+str(timeclim)+'Ma_'+CO2curve+'/sed.tif')

    envi_var = ['temp','salt','fstream','curr']
    for k in range(len(envi_var)):
        if envi_var[k] == 'temp' or envi_var[k] == 'salt':
            gtmp = shelfdata[envi_var[k]+'Min'].rio.write_nodata(np.nan, inplace=True).to_dataset()
            gtmp = xr.DataArray(
                    gtmp[envi_var[k]+'Min'].values,
                    dims=["latitude", "longitude"],
                    coords={"latitude": gtmp.lat.values[:,0], "longitude": gtmp.lon.values[0,:]},
                )
            gtmp = gtmp.to_dataset(name=envi_var[k]+'Min')
            gtmp.rio.write_crs('EPSG:4326', inplace=True)
            gtmp[envi_var[k]+'Min'].rio.to_raster('env_var/bio_res2_'+str(timeclim)+'Ma_'+CO2curve+'/'+envi_var[k]+"Min.tif") 

            gtmp = shelfdata[envi_var[k]+'Max'].rio.write_nodata(np.nan, inplace=True).to_dataset()
            gtmp = xr.DataArray(
                    gtmp[envi_var[k]+'Max'].values,
                    dims=["latitude", "longitude"],
                    coords={"latitude": gtmp.lat.values[:,0], "longitude": gtmp.lon.values[0,:]},
                )
            gtmp = gtmp.to_dataset(name=envi_var[k]+'Max')
            gtmp.rio.write_crs('EPSG:4326', inplace=True)
            gtmp[envi_var[k]+'Max'].rio.to_raster('env_var/bio_res2_'+str(timeclim)+'Ma_'+CO2curve+'/'+envi_var[k]+"Max.tif")
        else:
            gtmp = shelfdata[envi_var[k]].rio.write_nodata(np.nan, inplace=True).to_dataset()
            gtmp = xr.DataArray(
                    gtmp[envi_var[k]].values,
                    dims=["latitude", "longitude"],
                    coords={"latitude": gtmp.lat.values[:,0], "longitude": gtmp.lon.values[0,:]},
                )
            gtmp = gtmp.to_dataset(name=envi_var[k])
            gtmp.rio.write_crs('EPSG:4326', inplace=True)
            gtmp[envi_var[k]].rio.to_raster('env_var/bio_res2_'+str(timeclim)+'Ma_'+CO2curve+'/'+envi_var[k]+".tif")