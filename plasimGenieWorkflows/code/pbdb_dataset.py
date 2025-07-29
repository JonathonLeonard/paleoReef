import numpy as np
import pandas as pd
from code import paleoRotator

class PBDBDataset:
    def __init__(self, csv_path):
        """
        Initializes the PBDBDataset with the path to the dataset.

        :param data_path: Path to the PBDB dataset CSV file.
        """
        self.csv_path = csv_path
        self.df = self.load_and_clean()
        

    def load_and_clean(self):
        """
        Loads the PBDB dataset from the CSV file and performs basic cleaning.
        """
        self.df = pd.read_csv(self.csv_path, skiprows=list(range(0,18)))

        # keep the columns we need
        self.df = self.df[['occurrence_no', 'max_ma', 'min_ma', 'lng', 'lat','cc',
                            'paleomodel', 'geoplate', 'paleolng', 'paleolat',   
       'paleomodel2', 'geoplate2', 'paleolng2', 'paleolat2']]
        self.df['mid'] = round(((self.df['max_ma'] + self.df['min_ma']) / 2), 2)
        col = self.df.pop('mid')
        self.df.insert(0, 'mid', col)    
        self.df['rnd10'] = self.df['mid'].round(-1)
        # remove rows with duplicate lat/lng/age pairs
        self.df.drop_duplicates(subset=['lat', 'lng', 'rnd10'], inplace=True)
        # remove rows with missing lat/lng pairs
        self.df.dropna(subset=['lat', 'lng'], inplace=True)
        self.df = self.df.rename(columns={'lng':'longit', 'lat':'latit', 'paleomodel':'pmodel', 'paleolng':'plong',
                                          'paleolat':'plat', 'paleomodel2':'pmodel2', 'paleolng2':'plong2',
                                          'paleolat2':'plat2', 'geoplate':'gplate', 'geoplate2':'gplate2'})
        # replace where it says 'not computable using this model' with 'NaN'
        self.df.replace('not computable using this model', '??', inplace=True)
        # set the index to occurrence_no
        # self.df.set_index('occurrence_no', inplace=True)
        return self.df.reset_index(drop=True)
    
    def basic_df(self):
        df = self.df[['longit', 'latit', 'mid', 'rnd10', 'min_ma', 'max_ma', 'occurrence_no']]
        df['longit'] = df['longit'].round(2)
        df['latit'] = df['latit'].round(2)
        return df.rename(columns={'longit':'lon', 'latit':'lat', 'mid':'age', 'rnd10':'age10', 'min_ma':'minAge', 'max_ma':'maxAge', 'occurrence_no':'name'})
    
class paredDataset:
    """
    Placeholder for the paleoreef dataset processing
    """
    def __init__(self, csv_path):
        """
        Initializes the paredDataset with the path to the dataset.

        :param data_path: Path to the pared dataset CSV file.
        """
        self.csv_path = csv_path
        self.df = self.load_and_clean()

    def load_and_clean(self):
        """
        Loads the pared dataset from the CSV file and performs basic cleaning.
        """
        paredRaw = pd.read_csv(self.csv_path)
        # Do above for system y
        pared = paredRaw[paredRaw['top'] <= 255] # Remove any data older than what we're interested in
        pared = pared[(pared['biota_main'] == 1) | (pared['biota_sec']==1)] # These are the biota associated with warm water reefs
        pared = pared.drop(['systemCol', 'seriesCol','col'], axis=1) # Drop some random unnedded columns
        # remove rows where 'system.x' = 'Permian'
        pared = pared[pared['system.x'] != 'Permian'].reset_index(drop=True) # Remove any data older than what we're interested in

        pared['nearest10'] = pared['mid'].round(-1)

        self.df = pared
        return self.df.reset_index(drop=True)
    
    def basic_df(self):
        df = self.df[['longit', 'latit', 'mid', 'nearest10', 'min_ma', 'max_ma', 'r_number']]
        df['longit'] = df['longit'].round(2)
        df['latit'] = df['latit'].round(2)
        return df.rename(columns={'longit':'lon', 'latit':'lat', 'mid':'age', 'nearest10':'age10', 'min_ma':'minAge', 'max_ma':'maxAge', 'r_number':'name'})
    
class basicLocalityData:
    """
    Reconstruct and other things on consistently formatted locality data
    """
    def __init__(self, df):
        """
        Initializes the basicLocalityData with a DataFrame.

        :param df: DataFrame containing locality data.
        dataframe should have columns: 'lon', 'lat', 'age', 'age10', 'minAge', 'maxAge'
        """
        self.df = df

        # convert to 0-360 lons
        self.df['lon'] = self.df['lon'].apply(lambda x: x + 360 if x < 0 else x)
        
    def reconstruct_points(self, paleoRotatorObject, columnPrefix = '', anchor_plate_id=0):
        """
        Reconstructs the paleoreff points
        """
        maxAges = np.maximum(self.df['maxAge'], self.df['age10'] + 4.99)
        minAges = np.minimum(self.df['minAge'], self.df['age10'] - 5)
        gpml = paleoRotatorObject.points2gpml(
            self.df['lon'].values, 
            self.df['lat'].values, 
            maxAges=maxAges, 
            minAges=minAges, 
            names=self.df.index.values
        )
        # gpml.write('test2.gpml')
        self.df[f'{columnPrefix}PlateID'] = [int(feature.get_reconstruction_plate_id()) for feature in gpml]
        self.df[f'{columnPrefix}lats'], self.df[f'{columnPrefix}lons'] = paleoRotatorObject.rotatePoints(gpml, self.df['age'], anchor_plate_id=anchor_plate_id)
        # print(self.df)
        # For some reason the rotation puts coords back into the -180 to 180 range, so we need to convert them back to 0-360
        self.df[f'{columnPrefix}lons'] = self.df[f'{columnPrefix}lons'].apply(lambda x: x + 360 if x < 0 else x)
        return self.df
    
    def selectPoints4Ages(self, age):
        """
        Selects points for a given age.
        
        :param age: Age to select points for.
        :return: DataFrame with selected points.
        """
        selectedDF = self.df[(self.df['age'] > age -5) & (self.df['age'] < age + 5)]
        
        if selectedDF.empty:
            print(f"No points found for age {age}")
            return None
        return selectedDF.reset_index(drop=True)
    
    def selectPoints4Epochs(self, age):
        from pyrolite.util.time import Timescale
        ts = Timescale()

        epochs = ts.Epochs.copy().reset_index(drop=True)
        epochs['ageDiff'] = np.abs(epochs['MeanAge'] - age)
        closestEpoch = epochs.loc[epochs['ageDiff'].idxmin()]

        # print(closestEpoch)

        selectedDF = self.df[(self.df['age'] < closestEpoch['Start']) & (self.df['age'] > closestEpoch['End'])]

        if selectedDF.empty:
            print(f"No points found for age {age}")
            return None
        return selectedDF.reset_index(drop=True)
    
    def printClosestEpoch(self, age):
        from pyrolite.util.time import Timescale
        ts = Timescale()

        epochs = ts.Epochs.copy().reset_index(drop=True)
        epochs['ageDiff'] = np.abs(epochs['MeanAge'] - age)
        closestEpoch = epochs.loc[epochs['ageDiff'].idxmin()]

        print(f"Closest epoch to {age} Ma is {closestEpoch['Name']} ({closestEpoch['Start']} - {closestEpoch['End']} Ma) with a difference of {closestEpoch['ageDiff']} Ma")
        return closestEpoch

    @staticmethod
    def find_nearest_valid(ds, lon, lat, OGlon, OGlat, var, radius=5.6125):
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
            return np.nan, lon, lat
        # Find the coordinates of the valid point with the smallest distance
        min_dist_coords = dist.where(dist == dist.where(~np.isnan(valid_subset[var])).min(), drop=True)
        # print(min_dist_coords)
        if len(min_dist_coords.longitude) > 1:
            # choose the one closest to the original coordinates
            lons = min_dist_coords.longitude.values
            lons = np.abs(lons - OGlon)
            nearest_lon = min_dist_coords.longitude.values[np.where(lons == lons.min())].item()
        else:
            nearest_lon = min_dist_coords.longitude.values.item()
        nearest_lat = min_dist_coords.latitude.values.item()
        # Select the nearest valid value
        # print(nearest_lon, nearest_lat)
        nearest = valid_subset.sel(longitude=nearest_lon, latitude=nearest_lat, method='nearest')
        # return nearest[var].to_dataframe().reset_index()
        return nearest[var].values.item(), nearest_lon, nearest_lat
    
    @staticmethod
    def getNearestModelCellCoordinates(df, column_prefix, modelData):
        # modelData['landNAN'] = modelData['landsea'].where(modelData['landsea'] == 0, np.nan, modelData['landsea'])
        nearestValid = []
        df[f'{column_prefix}_nearestLons'] = np.nan
        df[f'{column_prefix}_nearestLats'] = np.nan
        for i, (lon,lat) in enumerate(zip(df[f'{column_prefix}_lons'], df[f'{column_prefix}_lats'])):
            # This will get the coords of the cell nearest to the reef point and whether the cell is land or ocean
            landseaVal = modelData['landNAN'].sel(longitude=lon, latitude=lat, method='nearest')
            # print(f'Finding nearest valid cell for {lon}, {lat} - landsea value: {landseaVal.values.item()}')
            if np.isnan(landseaVal.values.item()):
                # Find the nearest ocean cell - this also needs the non-interpolated reef coords as a 'tie breaker' when the nan cell is between two ocean cells
                nearest = basicLocalityData.find_nearest_valid(modelData, landseaVal.longitude.values.item(), landseaVal.latitude.values.item(), lon, lat, 'landNAN', radius=6)

            else:
                # nearest = (np.nan, np.nan, np.nan) # return nans if not within the defined radius
                nearest = (landseaVal.values.item(), landseaVal.longitude.values.item(), landseaVal.latitude.values.item())
            # print(f'Nearest valid cell for {lon}, {lat} is {nearest}')
            nearestValid.append(nearest)
            df.at[i, f'{column_prefix}_nearestLons'] = nearest[1]
            df.at[i, f'{column_prefix}_nearestLats'] = nearest[2]
        return df