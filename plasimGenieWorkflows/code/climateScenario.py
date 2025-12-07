import os
import importlib.util
import numpy as np
import pandas as pd
import pickle
import pygplates
import pygmt
import rioxarray as rio
import sys
import xarray as xr


class ClimateScenario:
    def __init__(self, climate_simulation_directory, CO2s, ages, topography_path=None):
        """
        Chooses the ideal climate scenario based on closestg match with 
        paleoenvironmental reconstructions.
        """
        self.climate_simulation_directory = climate_simulation_directory
        self.ages = ages
        self.CO2s = CO2s
        self.topoPath = topography_path

        self.paleotemperature_labels = []
        self.climate_simulations = {}
        for simulation in os.listdir(climate_simulation_directory):
            CO2 = simulation.split('_')[-1]
            if CO2 in self.CO2s:
                # self.climate_simulations.append(os.path.join(climate_simulation_directory, simulation))
                self.climate_simulations[CO2] = os.path.join(climate_simulation_directory, simulation)

        self.paleotemperature_df = pd.DataFrame(index=self.ages)
        self.SST_df = pd.DataFrame(index=self.ages)
        self.CO2_scenario_matches = pd.DataFrame(index=self.ages, columns=self.paleotemperature_labels)
        
        self.envRastersPath = f'data/biomodEnvRasters/{self.climate_simulation_directory.split("/")[-1]}'
        if not os.path.exists(self.envRastersPath):
            os.makedirs(self.envRastersPath)

    def process_merdith(self, file, pySCION_classes_path, ages=None, standalone=False, interpolate=True):
        """
        Process the Merdith et al. (2024) dataset to extract temperature values for specified ages.
        """
        if ages is None:
            ages = self.ages
        spec = importlib.util.spec_from_file_location("pySCION_classes", pySCION_classes_path)
        pySCION_classes = importlib.util.module_from_spec(spec)
        sys.modules["pySCION_classes"] = pySCION_classes
        spec.loader.exec_module(pySCION_classes)

        results = []
        with open(file, 'rb') as file_:
            results.append(pickle.load(file_))
        result = results[0]
        times = np.abs(np.mean(result[0].time_myr, axis=0))
        temps = np.mean(result[0].T_gast, axis=0)

        df = pd.DataFrame({'Age': times, 'paleoTemp': temps})
        if interpolate == False:
            return df[df['Age']<=255].set_index('Age').sort_index()
        agesDF = pd.DataFrame({'Age': ages}) # Create a dataframe with the ages we are interested in
        df = pd.concat([df, agesDF], ignore_index=True) # Concatenate the two dataframes
        df = df.set_index('Age').sort_index()

        df = df.interpolate(method='linear', axis=0).reset_index() # Interpolate the data to fill in missing values

        df = df[df['Age'].isin(ages)].set_index('Age')

        #remove duplicate first row if it exists and keep second row
        df = df[~df.index.duplicated(keep='last')]

        self.merdithLabel = 'MER25'
        self.paleotemperature_labels.append(self.merdithLabel)
        # join the dataframe with the paleotemperature_df

        df = df.rename(columns={'paleoTemp': self.merdithLabel})

        self.paleotemperature_df = self.paleotemperature_df.join(df, how='outer')

        if standalone:
            return df
        else:
            return self.paleotemperature_df

    def process_scotese(self, file, ages=None, standalone=False, tropical=False, interpolate=True):
        """
        Process the Scotese data to extract temperature values for specified ages.
        """
        if ages is None:
            ages = self.ages
        
        if interpolate == False:
            scoteseVar = pd.read_excel(file, '1my version')  # Read the specific sheet
        else:
            scoteseVar = pd.read_excel(file, '1my version')  # Read the specific sheet
        # df = scoteseVar[scoteseVar['Age'].isin(ages)][['Age', 'GAT']].set_index('Age').rename(columns={'GAT': 'paleoTemp'})
        if tropical:
            df = scoteseVar[['Age', 'Tropical']].rename(columns={'Tropical': 'paleoTemp'})
            self.scoteseLabel = "SCO21_TropicalAir"
        else:
            df = scoteseVar[['Age', 'GAT']].rename(columns={'GAT': 'paleoTemp'})
            self.scoteseLabel = "SCO21_GMST"

        if interpolate == False:
            return df[df['Age']<=255].set_index('Age').sort_index()
        agesDF = pd.DataFrame({'Age': ages}) # Create a dataframe with the ages we are interested in
        try:
            df = pd.concat([df, agesDF], ignore_index=True) # Concatenate the two dataframes
            df = df.drop_duplicates(subset='Age', keep='first').set_index('Age').sort_index()  # Remove duplicates, keeping the last occurrence

            df = df.interpolate(method='linear', axis=0).reset_index() # Interpolate the data to fill in missing values
            df = df[df['Age'].isin(ages)].set_index('Age')
        except Exception as e:
            print(f"Error processing Scotese data: {e}")

        
        self.paleotemperature_labels.append(self.scoteseLabel)
        df = df.rename(columns={'paleoTemp': self.scoteseLabel})
        self.paleotemperature_df = self.paleotemperature_df.join(df, how='outer')
        if standalone:
            return df
        else:
            return self.paleotemperature_df

    def process_phanDA(self, csv, ages=None, standalone=False, interpolate=True):
        """
        Process the PhanDA CSV file to extract temperature and CO2 data.
        """
        df = pd.read_csv(csv)
        if ages is None:
            ages = self.ages
        if interpolate == False:
            return df[df['AverageAge']<=255][['AverageAge', 'GMST_50']].set_index('AverageAge').sort_index().rename(columns={'GMST_50': 'paleoTemp'})
        agesDF = pd.DataFrame({'AverageAge': ages}) # Create a dataframe with the ages we are interested in
        df = pd.concat([df, agesDF], ignore_index=True) # Concatenate the two dataframes
        df = df.drop_duplicates(subset='AverageAge', keep='first').set_index('AverageAge').sort_index()  # Remove duplicates, keeping the last occurrence
        df = df.drop(columns=['Period','Epoch','Stage'])
        df = df.interpolate(method='linear', axis=0).reset_index() # Interpolate the data to fill in missing values

        df[df['AverageAge'] == 0] = df.iloc[1] # Make AverageAge = 0 the same values as the second row
        df.iloc[0,0] = 0 # Set the first age to 0

        # Just take the data we need
        df = df[df['AverageAge'].isin(ages)][['AverageAge', 'GMST_50']].set_index('AverageAge').rename(columns={'GMST_50': 'paleoTemp'})

        self.phanDAlabel = "PhanDA_GMST"
        self.paleotemperature_labels.append(self.phanDAlabel)
        df = df.rename(columns={'paleoTemp': self.phanDAlabel})
        self.paleotemperature_df = self.paleotemperature_df.join(df, how='outer')
        if standalone:
            return df
        else:
            return self.SST_df

    def process_phanDA_sst(self, excel, ages=None, standalone=False, interpolate=True):
        # For the tropical SST file from Emily Judd
        df = pd.read_excel(excel)
        if ages is None:
            ages = self.ages
        if interpolate == False:
            return df[df['AverageAge']<=255][['AverageAge', 'SST_30S30N_50']].set_index('AverageAge').sort_index().rename(columns={'SST_30S30N_50': 'paleoTemp'})
        agesDF = pd.DataFrame({'AverageAge': ages}) # Create a dataframe with the ages we are interested in
        df = pd.concat([df, agesDF], ignore_index=True) # Concatenate the two dataframes
        df = df.drop_duplicates(subset='AverageAge', keep='first').set_index('AverageAge').sort_index()  # Remove duplicates, keeping the last occurrence
        df = df.drop(columns=['Period','Epoch','Stage'])
        df = df.interpolate(method='linear', axis=0).reset_index() # Interpolate the data to fill in missing values
        df[df['AverageAge'] == 0] = df.iloc[1] # Make AverageAge = 0 the same values as the second row
        df.iloc[0,0] = 0 # Set the first age to 0
        # Just take the data we need
        df = df[df['AverageAge'].isin(ages)][['AverageAge', 'SST_30S30N_50']].set_index('AverageAge').rename(columns={'SST_30S30N_50': 'paleoTemp'})
        self.phanDAlabel = "PhanDA_SST"
        self.paleotemperature_labels.append(self.phanDAlabel)
        df = df.rename(columns={'paleoTemp': self.phanDAlabel})
        self.SST_df = self.SST_df.join(df, how='outer')

        if standalone:
            return df
        else:
            return self.SST_df

    def process_phanDA_sst_old(self, csv, ages=None, standalone=False):
        # For when we had just digitised the SST by eye
        df = pd.read_csv(csv)
        if ages is None:
            ages = self.ages
        agesDF = pd.DataFrame({'Age': ages}) # Create a dataframe with the ages we are interested in
        agesDF['Temp'] = np.nan
        df = pd.concat([df, agesDF], ignore_index=True).set_index('Age').sort_index() # Concatenate the two dataframes

        df = df.interpolate(method='linear', axis=0).reset_index() # Interpolate the data to fill in missing values

        # df[df['AverageAge'] == 0] = df.iloc[1] # Make AverageAge = 0 the same values as the second row
        # df.iloc[0,0] = 0 # Set the first age to 0

        # Just take the data we need
        df = df[df['Age'].isin(ages)][['Age', 'Temp']].set_index('Age').rename(columns={'Temp': 'paleoTemp'})

        self.phanDAlabel = "PhanDA_SST_paleotemps"
        self.paleotemperature_labels.append(self.phanDAlabel)
        df = df.rename(columns={'paleoTemp': self.phanDAlabel})
        self.SST_df = self.SST_df.join(df, how='outer')

        if standalone:
            return df
        else:
            return self.SST_df

    def process_songSST(self, file, ages=None, standalone=True, interpolate=True):
        if ages is None:
            ages = self.ages
        df = pd.read_excel(file, 'Table S1', skiprows=1)  # Read the specific sheet
        # # df = scoteseVar[scoteseVar['Age'].isin(ages)][['Age', 'GAT']].set_index('Age').rename(columns={'GAT': 'paleoTemp'})
        df = df[['Age (Ma)', 'Mean temperature (°C)']].rename(columns={'Mean temperature (°C)': 'SONG19'})
        if interpolate == False:
            return df[df['Age (Ma)']<=255].set_index('Age (Ma)').sort_index()
        agesDF = pd.DataFrame({'Age (Ma)': ages}) # Create a dataframe with the ages we are interested in
        df = pd.concat([df, agesDF], ignore_index=True) # Concatenate the two dataframes
        df = df.drop_duplicates(subset='Age (Ma)', keep='first').set_index('Age (Ma)').sort_index()  # Remove duplicates, keeping the last occurrence

        df = df.interpolate(method='linear', axis=0).reset_index() # Interpolate the data to fill in missing values
        df = df[df['Age (Ma)'].isin(ages)].set_index('Age (Ma)')

        self.songSSTlabel = "SONG19"
        self.paleotemperature_labels.append(self.songSSTlabel)
        df = df.rename(columns={'SONG19': self.songSSTlabel})
        
        return df

    def process_all_paleotemperatures(self, merdith_file=None, pySCION_classes_path=None, scotese_file=None, phanDA_file=None):
        """
        Run all paleotemperature processing methods and update self.paleotemperature_df.
        """
        if scotese_file:
            self.process_scotese(scotese_file)
        else:
            print("Scotese file not provided, skipping Scotese processing.")
        if merdith_file and pySCION_classes_path:
            self.process_merdith(merdith_file, pySCION_classes_path)
        else:
            print("Merdith file or pySCION classes path not provided, skipping Merdith processing.")
        if phanDA_file:
            self.process_phanDA(phanDA_file)
        else:
            print("PhanDA file not provided, skipping PhanDA processing.")

        if not merdith_file and not scotese_file and not phanDA_file:
            print("No paleotemperature data files provided, returning empty DataFrame.")

        return self.paleotemperature_df

    def process_grossman(self, file, ages=None, standalone=False, interpolate=True):
        """
        Process the Grossman et al. (2022) dataset for tropical SST proxies
        """
        if ages is None:
            ages = self.ages
        agesDF = pd.DataFrame({'Age [Ma]': ages})

        grossmanVar = pd.ExcelFile(file)
        grossmanVar = pd.read_excel(grossmanVar, 'Table S2', skiprows=1) # Read the specific sheet

        df = grossmanVar[['Age [Ma]', 'Locfit mean T (PO4+CaCO3; iceV, lat, arag, bel corr)']]
        if interpolate == False:
            return df[df['Age [Ma]']<=255].set_index('Age [Ma]').sort_index().rename(columns={'Locfit mean T (PO4+CaCO3; iceV, lat, arag, bel corr)': 'GRO22'})
        df = pd.concat([df, agesDF], ignore_index=True) # Concatenate the two dataframes

        df = df.drop_duplicates(subset='Age [Ma]', keep='first').set_index('Age [Ma]').sort_index()  # Remove duplicates, keeping the last occurrence
        df = df.interpolate(method='linear', limit_direction='both')
        df = df[df.index.isin(ages)][['Locfit mean T (PO4+CaCO3; iceV, lat, arag, bel corr)']].rename(columns={'Locfit mean T (PO4+CaCO3; iceV, lat, arag, bel corr)': 'GRO22'})

        self.SST_df = self.SST_df.join(df, how='outer')

        if standalone:
            return df
        else:
            return self.SST_df
        # return df

    def get_global_plasim_aves(self,var='puma_temperature_surface_air', tropical=False): 
        globalAveDF = pd.DataFrame(index=self.ages, columns=self.CO2s)
        for i, CO2 in enumerate(self.CO2s):
            model = self.climate_simulations[CO2]
            for age in self.ages:
                plasimGrid = f'{model}/{str(age)}Ma_nc_outputs/3000-12-30_plasim.nc'

                if not os.path.exists(plasimGrid):
                    print(f'{plasimGrid} does not exist, continuing but this may not work')
                    continue

                ds = xr.open_dataset(plasimGrid)
                airTemps = ds[var]

                if tropical:
                    airTemps = airTemps.sel(latitude=slice(18, -18))
                

                try:
                    info = pygmt.grdinfo(grid=airTemps, per_column=True, force_scan='2',verbose="e")
                    mean = round(float(info.split()[10]), 4)
                except Exception as e:
                    print(f"Error processing {plasimGrid}: {e}")
                    return airTemps

                globalAveDF.at[age,CO2] = mean 
        return globalAveDF
    
    def get_global_goldstein_aves(self,var='temp',tropical=False, song=False): 
        globalAveDF = pd.DataFrame(index=self.ages, columns=self.CO2s)
        for i, CO2 in enumerate(self.CO2s):
            model = self.climate_simulations[CO2]
            for age in self.ages:
                # print('\n')
                goldsteinGrid = f'{model}/{str(age)}Ma_nc_outputs/gold_spn_av_0000003000_00.nc'

                if not os.path.exists(goldsteinGrid):
                    print(f'{goldsteinGrid} does not exist, continuing but this may not work')
                    continue
                
                # Get the grid data
                ds = xr.open_dataset(goldsteinGrid)
                SSTs = ds[var][0][0]

                if tropical:
                    # Slice the grid to between lats of -30 and 30 degrees
                    SSTs = SSTs.sel(latitude=slice(-33, 33))
                if song:
                    # Slice the grid to between lats of -30 and 30 degrees
                    SSTs = SSTs.sel(latitude=slice(-43, 43))

                
                if song:
                    # For the songSST, we want the non-geographically weighted mean of the SSTs
                    mean = np.nanmean(SSTs.values)
                    # print(f'non-weighted mean: {mean}')
                    # info = pygmt.grdinfo(grid=SSTs, V="q", per_column=True, force_scan='2')
                    # mean = round(float(info.split()[10]), 4)
                    # print(f'mean from grdinfo: {mean}')
                else:
                    info = pygmt.grdinfo(grid=SSTs, per_column=True, force_scan='2',verbose="q")
                    mean = round(float(info.split()[10]), 4)

                globalAveDF.at[age,CO2] = mean 
        self.globalAvesGolstein = globalAveDF.copy()
        return self.globalAvesGolstein
    
    def get_paleotemperature_df(self, var='puma_surface_temperature_air', merdith_file=None, pySCION_classes_path=None, scotese_file=None, phanDA_file=None):
        self.globalAves = self.get_global_plasim_aves(var=var)
        self.process_all_paleotemperatures(merdith_file, pySCION_classes_path, scotese_file, phanDA_file)
        # return globalAves.join(ptempDF,how='outer')
        return self.globalAves, self.paleotemperature_df

    def get_matching_scenario(self, globalAvesDF, paleoVar_df=None):
        """
        Find the climate scenario that best matches the paleoenvironmental reconstructions.
        """
        if paleoVar_df is None or paleoVar_df.empty:
            print("Paleotemperature DataFrame is empty. Please run get_paleotemperature_df() first.")
            return None
        
        for paleotemp in paleoVar_df.columns:
            dfDiff = pd.DataFrame(index=self.ages)
            for CO2 in globalAvesDF.columns:
                if CO2 not in self.CO2s:
                    continue
                dfDiff[CO2] = np.abs(paleoVar_df[paleotemp] - globalAvesDF[CO2])

            for age in dfDiff.index:
                closestCO2 = dfDiff.loc[age].idxmin()
                self.CO2_scenario_matches.at[age, paleotemp] = closestCO2

        return self.CO2_scenario_matches

    def get_environmental_rasters(self, for_fuzzy_logic=False):
        """
        Get the environmental rasters for use in biomod SDMs.
        """
        for i, CO2 in enumerate(self.climate_simulations):
            model = self.climate_simulations[CO2]
            for age in self.ages:
                plasimGrid = f'{model}/{str(age)}Ma_nc_outputs/3000-12-30_plasim.nc'
                goldsteinGrid = f'{model}/{str(age)}Ma_nc_outputs/gold_spn_av_0000003000_00.nc'

                try:
                    with xr.open_dataset(goldsteinGrid) as oceanDS:
                        atmosDS = xr.open_dataset(plasimGrid).interp_like(oceanDS) # Re-cast the plasim coordinates into the same as goldstein
                        # Variables into xarray
                        # Goldstein
                        salinity = oceanDS['salinity'][0][0]
                        salinity = xr.where(salinity < 20, np.nan, salinity) # To remove negative/very low values
                        upwelling = oceanDS['wvel'][0][1]
                        # solarFlux = oceanDS['netsolar'][0]
                        runoff = oceanDS['runoff'][0]
                        bathymetry = oceanDS['bathymetry']

                        # Plasim
                        landsea = atmosDS['grid_landsea_mask'] 
                        landseaNAN = xr.where(landsea > 0.1, np.nan, 1) # Impose NaNs on land
                        surfaceTemp = atmosDS['puma_temperature_surface'] * landseaNAN
                        SST_DJF = (atmosDS['DJF_GENIE_SST'] - 273.15) * landseaNAN
                        SST_JJA = (atmosDS['JJA_GENIE_SST'] - 273.15) * landseaNAN
                        solar = atmosDS['ebal_surface_solar_radiation'] * landseaNAN
                        solar_DJF = atmosDS['DJF_incoming_solar'] * landseaNAN
                        solar_JJA = atmosDS['JJA_incoming_solar'] * landseaNAN

                        # Combine seasonal data to get summer and winter means
                        # SST
                        SST_DJF_northern = SST_DJF.where(SST_DJF.latitude < 0, 0)
                        SST_DJF_southern = SST_DJF.where(SST_DJF.latitude > 0, 0)
                        SST_JJA_northern = SST_JJA.where(SST_JJA.latitude < 0, 0)
                        SST_JJA_southern = SST_JJA.where(SST_JJA.latitude > 0, 0)
                        SST_summer = SST_DJF_northern + SST_JJA_southern
                        SST_winter = SST_DJF_southern + SST_JJA_northern
                        # Solar Insolation
                        solar_DJF_northern = solar_DJF.where(solar_DJF.latitude < 0, 0)
                        solar_DJF_southern = solar_DJF.where(solar_DJF.latitude > 0, 0)
                        solar_JJA_northern = solar_JJA.where(solar_JJA.latitude < 0, 0)
                        solar_JJA_southern = solar_JJA.where(solar_JJA.latitude > 0, 0)
                        solar_summer = solar_DJF_northern + solar_JJA_southern
                        solar_winter = solar_DJF_southern + solar_JJA_northern
                        
                        # Remove white spaces from variable names
                        # netsolar.attrs['long_name'] = 'netSolarHeatFlux'
                        bathymetry.attrs['long_name'] = 'oceanDepth'
                        upwelling.attrs['long_name'] = 'verticalCurrent'

                        if not for_fuzzy_logic:
                            # Save as a geotiff
                            save_path = f'{self.envRastersPath}/{CO2}ppm/{age}Ma_envRaster'
                            if not os.path.exists(save_path):
                                os.makedirs(save_path)
                            surfaceTemp.rio.to_raster(f'{save_path}/SST_annual.tif')
                            SST_summer.rio.to_raster(f'{save_path}/SST_summer.tif')
                            SST_winter.rio.to_raster(f'{save_path}/SST_winter.tif')
                            solar.rio.to_raster(f'{save_path}/solar_annual.tif')
                            solar_summer.rio.to_raster(f'{save_path}/solar_summer.tif')
                            solar_winter.rio.to_raster(f'{save_path}/solar_winter.tif')
                            salinity.rio.to_raster(f'{save_path}/salinity_annual.tif')
                            upwelling.rio.to_raster(f'{save_path}/upwelling_annual.tif')
                            runoff.rio.to_raster(f'{save_path}/runoff_annual.tif')
                            bathymetry.rio.to_raster(f'{save_path}/bathymetry.tif')
                        else:
                            save_path = f'{self.envRastersPath}/SSTs/{CO2}ppm/{age}Ma'
                            if not os.path.exists(save_path):
                                os.makedirs(save_path)
                            atmosDS['SST_summer'] = SST_summer
                            atmosDS['SST_winter'] = SST_winter
                            # save as netcdf
                            # Remove the FillValue that xarray automaticcaly puts in (stuffs up some of the plotting)
                            enc = {'_FillValue': None}
                            encoding = {var:enc for var in atmosDS.variables}
                            atmosDS.to_netcdf(f'{save_path}/plasim_outs_SSTs.nc', encoding=encoding)
                        oceanDS.close()
                        atmosDS.close()
                except Exception as e:
                    print(f"Error processing {plasimGrid} or {goldsteinGrid}: {e}")
                    continue

    # @staticmethod            
    # def getSSTProbs(ds, varName = 'temp', tempAbsMin=18, tempMin=23, tempMax=28, tempAbsMax=32):
    #     """
    #     Get the fuzzy logic-like probabilities of SST with exponential decay from boundaries
    #     """
    #     # Assign all temp values between 23 and 28 °C to 1
    #     ds['tempProb'] = xr.where(
    #         (ds[varName] > tempMin) & (ds[varName] < tempMax),
    #         1,
    #         0
    #     )
    #     # Give an exponentially declining value from 29 to 31 °C to 0
    #     ds['tempProb'] = xr.where(
    #         (ds[varName] > tempMax) & (ds[varName] < tempAbsMax),
    #         np.exp(-((ds[varName] - tempMax) / 2)),
    #         ds['tempProb']
    #     )

    #     # Same for 20 to 23 °C
    #     ds['tempProb'] = xr.where(
    #         (ds[varName] > tempAbsMin) & (ds[varName] < tempMin),
    #         np.exp((ds[varName] - tempMin) / 2),
    #         ds['tempProb']
    #     )
        
    #     return ds

    @staticmethod
    def getSSTProbs(ds, varName='temp',tempMin=22, tempMax=30):
        """
        Get the fuzzy logic-like probabilities of SST with exponential decay from boundaries
        """
        # Assign all temp values between 23 and 28 °C to 1
        ds['tempProb'] = xr.where(
            (ds['SST_winter'] > tempMin) & (ds['SST_summer'] < tempMax),
            1,
            0
        )

        ds['tempProb'] = xr.where(
            ds['SST_summer'] > tempMax,
            np.exp(-(ds['SST_summer'] - tempMax) / 3),
            ds['tempProb']
        )

        ds['tempProb'] = xr.where(
            ds['SST_winter'] < tempMin,
            np.exp((ds['SST_winter'] - tempMin) / 5),
            ds['tempProb']
        )
        
        return ds

    def fuzzy_logic_projection(self, matching_scenario=None, CO2=None):
        """
        Run the fuzzy logic projection using the specified model.
        """
        if matching_scenario is None and CO2 is None:
            print("No matching scenario or CO2 provided, returning empty dictionary.")
            return {}
        if matching_scenario is not None and CO2 is not None:
            print("Both matching scenario and CO2 provided, using matching scenario.")
            CO2 = None

        climateSims = {}
        for age in self.ages:
            if CO2 is not None:
                co2 = CO2
            else:
                co2 = matching_scenario[age]

            for folder in os.listdir(self.climate_simulation_directory):
                if folder.split('_')[-1] == co2:
                    climateData = f'{self.climate_simulation_directory}/{folder}/{age}Ma_nc_outputs/gold_spn_av_0000003000_00.nc'
                    # climateData = f'{self.climate_simulation_directory}/{folder}/{age}Ma_nc_outputs/3000-12-30_plasim.nc'
                    try:
                        ds = xr.open_dataset(climateData)
                        SST_path = xr.open_dataset(f'{self.envRastersPath}/SSTs/{co2}ppm/{age}Ma/plasim_outs_SSTs.nc')
                    except FileNotFoundError as e:
                        print(f"File not found: {e}")
                        continue

                    SST_summer = SST_path['SST_summer']
                    SST_winter = SST_path['SST_winter']


                    ds['SST_summer'] = SST_summer
                    ds['SST_winter'] = SST_winter
                    climateSims[age] = ds

        SSTprobDatasets = {}
        for age, ds in climateSims.items():
            try:
               
                ds = self.getSSTProbs(ds, varName='temp')
                # ds.to_netcdf(f'{self.envRastersPath}/{matching_scenario[age]}ppm/{age}Ma_envRaster/SST_probabilities.nc')
            except Exception as e:
                print(f"Error processing {ds}: {e}")
                continue
            SSTprobDatasets[age] = ds
        return SSTprobDatasets
    
    



        

