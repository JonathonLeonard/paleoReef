import cv2 
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import pygmt    
import xarray as xr
import rioxarray as rio 
from scipy import interpolate
from scipy.interpolate import NearestNDInterpolator

class present_day_photozoan_presences:
    def __init__(self, presences_file):
        self.presences_file = presences_file

        # The size of the carbonate factory image in pixels
        self.carbfactory = np.zeros((797, 1595))

        self.elevationGrid = pygmt.datasets.load_earth_relief(resolution="10m",registration="gridline")

        self.img = cv2.imread(self.presences_file) 

    def load_presences(self):
        """
        Load the present-day photozoan presences from the Laugie 2019 image
        """
        self.img = cv2.imread(self.presences_file) 

    def mask_carbonates(self, photozoan_t_rgb_range = ((0,0,120), (70,65,200)), heterozoan_c_rgb_range = ((130,50,0), (190,135,90)), photozoan_c_rgb_range = ((85,0,55), (145,75,130)), biochemical_rgb_range = ((0,130,180), (80,170,240))):
        """
        Mask the carbonates in the image based on the RGB ranges provided.
        
        :param img: The image to mask.
        :param photozoan_t_rgb_range: RGB range for photozoan T carbonates.
        :param heterozoan_c_rgb_range: RGB range for heterozoan C carbonates.
        :param photozoan_c_rgb_range: RGB range for photozoan C carbonates.
        :param biochemican_rgb_range: RGB range for biochemical carbonates.
        """
        mask_photozoan_t = cv2.inRange(self.img, photozoan_t_rgb_range[0], photozoan_t_rgb_range[1])
        res_photozoan_t = cv2.bitwise_and(self.img, self.img, mask= mask_photozoan_t)

        mask_heterozoan_c = cv2.inRange(self.img, heterozoan_c_rgb_range[0], heterozoan_c_rgb_range[1])
        res_heterozoan_c = cv2.bitwise_and(self.img, self.img, mask= mask_heterozoan_c)

        mask_photozoan_c = cv2.inRange(self.img, photozoan_c_rgb_range[0], photozoan_c_rgb_range[1])
        res_photozoan_c = cv2.bitwise_and(self.img, self.img, mask= mask_photozoan_c)

        mask_biochem = cv2.inRange(self.img, biochemical_rgb_range[0], biochemical_rgb_range[1])
        res_biochem = cv2.bitwise_and(self.img, self.img, mask= mask_biochem)



        # Assign unique IDs to each type so we can combine inot one figure
        ids = np.where(np.logical_or(res_heterozoan_c[:,:,2]>0,res_heterozoan_c[:,:,1]>0))
        self.carbfactory[ids] = 4

        ids = np.where(np.logical_or(res_photozoan_c[:,:,2]>0,res_photozoan_c[:,:,1]>0))
        self.carbfactory[ids] = 3

        ids = np.where(np.logical_or(res_biochem[:,:,2]>0,res_biochem[:,:,1]>0))
        self.carbfactory[ids] = 1

        ids = np.where(np.logical_or(res_photozoan_t[:,:,2]>0,res_photozoan_t[:,:,1]>0))
        self.carbfactory[ids] = 2

        # Widen the figure
        ncarbfactory = np.zeros((801, 1595))

        ncarbfactory[2:-2,:] = self.carbfactory

        # Impose onto a geographical grid
        lon = np.linspace(-180,180,1595)
        lat = np.linspace(-90,90,801)

        mlon, mlat = np.meshgrid(lon, lat)
        interp = NearestNDInterpolator(list(zip(mlon.flatten(), mlat.flatten())), ncarbfactory.flatten())
        mlon2, mlat2 = np.meshgrid(self.elevationGrid.lon.values, self.elevationGrid.lat.values)
        self.Carbfact = interp(mlon2, mlat2)

        return self.Carbfact
    
    def assign_to_xarray(self):
        """
        Assign the carbonate factory data to an xarray DataArray.
        """
        ds_carbfact = xr.Dataset({})

        ds_carbfact['elevation'] = self.elevationGrid # This is the elevation grid downloaded at the start
        ds_carbfact['carbfac'] = (('lat','lon'),np.flipud(self.Carbfact))

        # Grid of carbonates where everything else is a NaN
        factory = ds_carbfact.carbfac.where(ds_carbfact.elevation<2000) # Also remove carbonates at high elevations
        factory2 = factory.where(factory>0.1)

        dsCarb = xr.Dataset({})
        dsCarb['elevation'] = self.elevationGrid
        dsCarb['carbfac'] = (('lat','lon'),factory2.values)

        dsCarb['photozoan'] = (('lat','lon'),dsCarb.carbfac.where(dsCarb.carbfac==2).values)
        dsCarb['photoC'] = (('lat','lon'),dsCarb.carbfac.where(dsCarb.carbfac==3).values)
        dsCarb['biochemical'] = (('lat','lon'),dsCarb.carbfac.where(dsCarb.carbfac==1).values)
        dsCarb['heterozoan'] = (('lat','lon'),dsCarb.carbfac.where(dsCarb.carbfac==4).values)

        self.dsCarb = dsCarb
        return self.dsCarb
    
    def save_photozoan_to_csv(self, output_dir='data/photozoan_locations'):
        photozoan = self.dsCarb.photozoan.where(self.dsCarb.photozoan != 2, 1)
        photozoanDF = photozoan.to_dataframe().dropna(axis=0).reset_index()
        photozoanDF.to_csv(f'{output_dir}/photozoan.csv') 

        self.photozoanDF = photozoanDF
        return self.photozoanDF

    def save_filtered_photozoan_to_csv(self, presentday_goldstein_file, output_dir='data/photozoan_locations'):
        """
        Filter to photozoans in shallow marine settings in the climate model
        """
        oceanDS = xr.open_dataset(presentday_goldstein_file)
        bathymetry = oceanDS['bathymetry']

        reefPts = self.photozoanDF.dropna(axis=0).reset_index(drop=True)

        reefPts = reefPts[(reefPts['lat']>-50) & (reefPts['lat']<50)] # Filter out points well outside the expected range
        # print(reefPts.head())
        elevationReefs = pygmt.grdtrack(bathymetry,reefPts, newcolname='bathymetry', incols=[1,0])
        elevationReefs = elevationReefs.rename(columns={'lat':'lon','lon':'lat'})
        # re-order the columns
        elevationReefs = elevationReefs[['lat','lon','photozoan']]

        # print(elevationReefs.head())

        filteredReefs = elevationReefs.dropna(axis=0)
        filteredReefs = filteredReefs[filteredReefs['photozoan']<1500]
        filteredReefs['photozoan'] = 1

        self.filteredReefs = filteredReefs

        filteredReefs.to_csv(f'{output_dir}/photozoan_filtered.csv',index=False)

        return self.filteredReefs

    def find_nearest_valid(self, ds, lon, lat, OGlon, OGlat, var, radius):
        """Find the nearest valid value in a dataset for a given variable and coordinates.
        Taken from https://github.com/pydata/xarray/issues/644"""
        # Select a subset of points within a certain radius
        subset = ds.sel(longitude=slice(lon-radius, lon+radius), latitude=slice(lat-radius, lat+radius))
        # Calculate the Euclidean distance to the points in the subset
        # dist = np.sqrt((subset.longitude - lon)**2 + (subset.latitude - lat)**2)
        dist = np.sqrt((subset.longitude - OGlon)**2 + (subset.latitude - OGlat)**2)
        # subset['landsea'].plot()

        # Create a new dataset that only includes valid values
        try:
            valid_subset = subset.where(~np.isnan(subset[var]), drop=True)
        except ValueError:
            print(f'No valid subset found at {lon}, {lat}')
            return None
        # Find the coordinates of the valid point with the smallest distance
        min_dist_coords = dist.where(dist == dist.where(~np.isnan(valid_subset[var])).min(), drop=True)
        # print(min_dist_coords)
        if len(min_dist_coords.longitude) > 1:
            # choose the one closest to the original coordinates
            lons = min_dist_coords.longitude.values
            lons = np.abs(lons - OGlon)
            if lons[0]-lons[1] < 0.01:
                # if the two closest points are very close together, just choose the first one
                lons = lons[0]
            nearest_lon = min_dist_coords.longitude.values[np.where(lons == lons.min())].item()

        else:
            nearest_lon = min_dist_coords.longitude.values.item()
        if len(min_dist_coords.latitude) > 1:
            # choose the one closest to the original coordinates
            lats = min_dist_coords.latitude.values
            lats = np.abs(lats - OGlat)
            nearest_lat = min_dist_coords.latitude.values[np.where(lats == lats.min())[0][0]].item()
        else:   
            nearest_lat = min_dist_coords.latitude.values.item()
        # Select the nearest valid value
        # print(nearest_lon, nearest_lat)
        nearest = valid_subset.sel(longitude=nearest_lon, latitude=nearest_lat, method='nearest')
        # return nearest[var].to_dataframe().reset_index()
        return nearest[var].values.item(), nearest_lon, nearest_lat
    
    def save_photozoan_to_csv_with_nearest(self, presentday_goldstein_file, output_dir='data/photozoan_locations', filtered_shallow_reefs_only = False):
        """
        Filter to photozoans in shallow marine settings in the climate model
        and find the nearest valid value for each point.
        """
        climateData = xr.open_dataset(presentday_goldstein_file)
        climateData['landsea'] = xr.where(climateData['bathymetry'] > 0, 0., 1)
        climateData['landNAN'] = climateData['landsea'].where(climateData['landsea'] == 0) # near variable where land values are NaN (instead of 1)

        if filtered_shallow_reefs_only:
            select = pd.DataFrame({
                'longitude': self.filteredReefs['lon'],
                'latitude': self.filteredReefs['lat']
            })
        else:
            select = pd.DataFrame({
            # 'longitude': np.where(df[paleoLon] >= 0, df[paleoLon], df[paleoLon] + 360),
            'longitude': self.photozoanDF['lon'],
            'latitude': self.photozoanDF['lat']
            })

        select['longitude'] = select['longitude'].where(select['longitude'] >= 0, select['longitude'] + 360) # convert to 0-360
        # drop dunplicate points
        select = select.drop_duplicates(subset=['longitude', 'latitude']).reset_index(drop=True)

        # Get cell coordinates of each reef point then move any points on land to the nearest ocean cell
        reefLandsea = []
        nearestValid = []
        for index, row in select.iterrows():
            # This will get the coords of the cell nearest to the reef point and whether the cell is land or ocean
            landseaVal = climateData['landNAN'].sel(longitude=row['longitude'], latitude=row['latitude'], method='nearest')
            reefLandsea.append((landseaVal.values.item(), landseaVal.longitude.values.item(), landseaVal.latitude.values.item()))

            if np.isnan(landseaVal.values.item()):
                # Find the nearest ocean cell - this also needs the non-interpolated reef coords as a 'tie breaker' when the nan cell is between two ocean cells
                nearest = self.find_nearest_valid(ds=climateData, lon=landseaVal.longitude.values.item(), lat=landseaVal.latitude.values.item(), OGlon=row['longitude'], OGlat=row['latitude'], var='landNAN', radius=6)
            else:
                nearest = (np.nan, np.nan, np.nan) # return nans if not within the defined radius
            nearestValid.append(nearest)

        # Now put the nearest coords and land/ocean values into the select dataframe
        select['reefLandsea'] = [landseaVal[0] if landseaVal is not None else np.nan for landseaVal in reefLandsea]
        select['realignedLon'] = [landseaVal[1] if landseaVal is not None else np.nan for landseaVal in reefLandsea]
        select['realignedLat'] = [landseaVal[2] if landseaVal is not None else np.nan for landseaVal in reefLandsea]
        select['nearestLandsea'] = [item[0] if item is not None else np.nan for item in nearestValid]
        select['nearestLon'] = [item[1] if item is not None else np.nan for item in nearestValid]
        select['nearestLat'] = [item[2] if item is not None else np.nan for item in nearestValid]

        # Have a final set of columns that is either the cell coord that the reef point is in or the one that is nearest to it
        select = select[(select['nearestLandsea'] == 0) | (select['reefLandsea'] == 0)]
        select['finalLon'] = np.where(select['nearestLandsea'] == 0, select['nearestLon'], select['realignedLon'])
        select['finalLat'] = np.where(select['nearestLandsea'] == 0, select['nearestLat'], select['realignedLat'])

        # Save out the coords of the centroid of present cells and save in same format as the other photozoan datasets
        onlyCellCoordsDF = select[['finalLat', 'finalLon']].drop_duplicates().reset_index(drop=True)
        onlyCellCoordsDF = onlyCellCoordsDF.rename(columns={'finalLat': 'lat', 'finalLon': 'lon'})
        onlyCellCoordsDF['photozoan'] = 1
        onlyCellCoordsDF.to_csv(f'{output_dir}/photozoan_nearest_modelCoords.csv', index=False)

        return select