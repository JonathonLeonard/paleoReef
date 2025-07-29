import os
import importlib.util
import numpy as np
import pandas as pd
import pickle
import pygplates
import pygmt
import shutil
import subprocess
import sys

class BiomodSDM:
    def __init__(self, env_rasters_path, training_data_path, SDM_data_path, respName, matching_scenario):
        self.envRastersPath = env_rasters_path
        self.modelName = self.envRastersPath.split('/')[-1]
        self.trainingDataPath = training_data_path
        self.SDMDataPath = SDM_data_path
        self.respName = respName
        self.matchingScenario = matching_scenario
        self.envRastersPaths = []
        self.presentDayRastersList = ''
        self.envRastersList = {}

    def load_env_raster_paths(self, matching_scenario_dict):
        """
        Load environmental rasters based on the matching scenario dictionary.
        
        :param matching_scenario_dict: Dictionary with age as keys and CO2 levels as values.
        :return: List of raster paths.
        """
        for age, CO2 in matching_scenario_dict.items():
            raster_path = f'{self.envRastersPath}/{CO2}ppm/{age}Ma_envRaster'
            if os.path.exists(raster_path):
                self.envRastersPaths.append(raster_path)
            else:
                print(f'Raster path {raster_path} does not exist, skipping.')
        
        return self.envRastersPaths
    
    def load_present_env_rasters(self, variables_list):
        """
        Load present-day environmental rasters based on the provided variables list.
        
        :param variables_list: List of environmental variable names.
        :return: List of raster paths.
        """
        for var in variables_list:
            raster_path = f'{self.envRastersPath}/280ppm/0Ma_envRaster/{var}.tif'
            if os.path.exists(raster_path):
                # self.presentDayRastersList.append(raster_path)
                self.presentDayRastersList = raster_path + ',' + self.presentDayRastersList
            else:
                print(f'Raster path {raster_path} does not exist, skipping.')
        
        return self.presentDayRastersList
    
    def load_past_env_rasters_list(self, variables_list, ages):
        """
        Load past environmental rasters based on the provided variables list and ages.
        
        :param variables_list: List of environmental variable names.
        :param ages: List of ages to consider.
        :return: List of raster paths.
        """
        for age in ages:
            tmp_rasters = ''
            for var in variables_list:
                raster_path = f'{self.envRastersPath}/{self.matchingScenario[age]}ppm/{age}Ma_envRaster/{var}.tif'
                if os.path.exists(raster_path):
                    tmp_rasters = raster_path + ',' + tmp_rasters
                else:
                    print(f'Raster path {raster_path} does not exist, skipping.')
            self.envRastersList[age] = tmp_rasters
        
        return self.envRastersList
        
    def run_biomod_ensemble_model(self):
        if self.presentDayRastersList == []:
            print('No environmental rasters found, exiting.')
            sys.exit(1)
        subprocess.call(['Rscript','code/getEnsembleModel.r', self.trainingDataPath, self.respName, self.presentDayRastersList], shell=False)
        # shutil.move(self.respName, self.SDMDataPath)
    
    def run_biomod_ensemble_projection(self, ages):
        for age in ages:
            co2 = self.matchingScenario[age]
            envRasters = self.envRastersList[age]
            subprocess.call(['Rscript','code/getEnsembleProjection.r', envRasters, str(age), co2, self.modelName], shell=False)


