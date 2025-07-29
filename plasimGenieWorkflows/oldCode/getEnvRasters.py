#!/usr/bin/python3

import matplotlib.pyplot as plt
import netCDF4 as nc
import numpy as np
import os
import pandas as pd
import pygmt
import rioxarray as rio
import shutil
import xarray as xr
import sys

def getEnvRasters(age,dpath):
	ocean = f'{dpath}/{age}Ma_nc_outputs/gold_spn_av_0000003000_00.nc'
	atmos = f'{dpath}/{age}Ma_nc_outputs/3000-12-30_plasim.nc'


	try:	
		with xr.open_dataset(ocean) as oceanDS:
			# oceanDS = xr.open_dataset(ocean).interp_like(atmosDS) #Re-cast the goldstein coordinates into the same as plasim
			atmosDS = xr.open_dataset(atmos).interp_like(oceanDS)
			# Variables into xarray
			# Goldstein
			salinity = oceanDS['salinity'][0][0]
			salinity = xr.where(salinity<20,np.nan,salinity) # To remove negative values
			upwelling = oceanDS['wvel'][0][1]
			# solarFlux = oceanDS['netsolar'][0]
			runoff = oceanDS['runoff'][0]
			bathymetry = oceanDS['bathymetry']

			# Plasim
			landsea = atmosDS['grid_landsea_mask']
			landseaNAN = xr.where(landsea>0.1, np.nan, 1) # To impose NaNs on land

			surfaceTemp = atmosDS['puma_temperature_surface']*landseaNAN
			SST_DJF = (atmosDS['DJF_GENIE_SST'] - 273.15)*landseaNAN
			SST_JJA = (atmosDS['JJA_GENIE_SST'] - 273.15)*landseaNAN
			solar = atmosDS['ebal_surface_solar_radiation']*landseaNAN
			solar_DJF = atmosDS['DJF_incoming_solar']*landseaNAN
			solar_JJA = atmosDS['JJA_incoming_solar']*landseaNAN

			oceanDS.close()
			atmosDS.close()

			# Combine seasonal data to get summer and winter means
			# SST
			SST_DJF_northern = SST_DJF.where(SST_DJF.latitude<0,0)
			SST_DJF_southern = SST_DJF.where(SST_DJF.latitude>0,0)

			SST_JJA_northern = SST_JJA.where(SST_JJA.latitude<0,0)
			SST_JJA_southern = SST_JJA.where(SST_JJA.latitude>0,0)

			SST_summer = SST_DJF_northern + SST_JJA_southern
			SST_winter = SST_DJF_southern + SST_JJA_northern

			# Solar Insolation
			solar_DJF_northern = solar_DJF.where(solar_DJF.latitude<0,0)
			solar_DJF_southern = solar_DJF.where(solar_DJF.latitude>0,0)

			solar_JJA_northern = solar_JJA.where(solar_JJA.latitude<0,0)
			solar_JJA_southern = solar_JJA.where(solar_JJA.latitude>0,0)

			solar_summer = solar_DJF_northern + solar_JJA_southern
			solar_winter = solar_DJF_southern + solar_JJA_northern

			# Remove white space from meta data
			# netsolar.attrs['long_name'] = 'netSolarHeatFlux'
			bathymetry.attrs['long_name'] = 'oceanDepth'
			upwelling.attrs['long_name'] = 'verticalCurrent'

			# Convert to geoTiff
			dpathStub = dpath.split('/')[-1] # Should be the model name
			envRasterPath = f'data/rasters/{dpathStub}/{age}Ma_envRaster'
			if os.path.exists(envRasterPath):
			    shutil.rmtree(envRasterPath)
			os.makedirs(envRasterPath)

			surfaceTemp.rio.to_raster(f'{envRasterPath}/SST_annual.tif')
			SST_summer.rio.to_raster(f'{envRasterPath}/SST_summer.tif')
			SST_winter.rio.to_raster(f'{envRasterPath}/SST_winter.tif')

			solar.rio.to_raster(f'{envRasterPath}/solar_annual.tif')
			solar_summer.rio.to_raster(f'{envRasterPath}/solar_summer.tif')
			solar_winter.rio.to_raster(f'{envRasterPath}/solar_winter.tif')

			# bathymetry.rio.to_raster(f'{envRasterPath}/bathymetry.tif')
			salinity.rio.to_raster(f'{envRasterPath}/salinity_annual.tif')
			upwelling.rio.to_raster(f'{envRasterPath}/upwelling_annual.tif')
			runoff.rio.to_raster(f'{envRasterPath}/runoff_annual.tif')
	except Exception as e:
		print(e)
		return None

def main():
	ages = np.arange(0,260,10)
	co2s=['100','200','280','560','1120','2240','4480','8960','17920','fosterCO2']
	os.chdir('/Users/jono/OneDrive - The University of Sydney (Staff)/paleoReef/plasimGenieWorkflows')

	if len(sys.argv)>1:
		dpathStub=sys.argv[1]
	else:
		# dpathStub='2024-10-08_realigned_scotese_paleo_shallowSea'
		dpathStub='2024-07-23_realigned_earthbyte_paleo_PaleomagRef'

	for co2 in co2s:
		dpath=f'../climateModelData/{dpathStub}_{co2}'
		# dpath=f'/Users/jono/OneDrive - The University of Sydney (Staff)/paleoReef/climateModelData/{dpathStub}_{co2}'
		for age in ages:
			print(age)
			getEnvRasters(age,dpath)
			print('\n')

if __name__ == "__main__":
	main()