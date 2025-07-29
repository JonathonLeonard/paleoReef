import importlib.util
import pandas as pd
import numpy as np
import pickle
import sys

# Some functions to load and clean different paleotemperature datasets

def process_merdith(file, ages, pySCION_classes_path):
    """
    Process the Merdith et al. (2024) dataset to extract temperature values for specified ages.
    """
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
    agesDF = pd.DataFrame({'Age': ages}) # Create a dataframe with the ages we are interested in
    df = pd.concat([df, agesDF], ignore_index=True) # Concatenate the two dataframes
    df = df.set_index('Age').sort_index()

    df = df.interpolate(method='linear', axis=0).reset_index() # Interpolate the data to fill in missing values

    df = df[df['Age'].isin(ages)].set_index('Age')

    #remove duplicate first row if it exists and keep second row
    df = df[~df.index.duplicated(keep='last')]

    label = 'Merdith2024_pyscion'
    return df, label

def process_scotese(file, ages):
    """
    Process the Scotese data to extract temperature values for specified ages.
    """
    scoteseVar = pd.read_excel(file, '1my version')  # Read the specific sheet
    df = scoteseVar[scoteseVar['Age'].isin(ages)][['Age', 'GAT']].set_index('Age').rename(columns={'GAT': 'paleoTemp'})
    
    label = "Scotese2021_paleotemps"
    return df, label

def process_phanDA(csv, ages):
    """
    Process the PhanDA CSV file to extract temperature and CO2 data.
    """
    df = pd.read_csv(csv)
    agesDF = pd.DataFrame({'AverageAge': ages}) # Create a dataframe with the ages we are interested in
    df = pd.concat([df, agesDF], ignore_index=True) # Concatenate the two dataframes
    df = df.set_index('AverageAge').sort_index() # Set the index to the age column
    df = df.drop(columns=['Period','Epoch','Stage'])
    df = df.interpolate(method='linear', axis=0).reset_index() # Interpolate the data to fill in missing values

    df[df['AverageAge'] == 0] = df.iloc[1] # Make AverageAge = 0 the same values as the second row
    df.iloc[0,0] = 0 # Set the first age to 0

    # Just take the data we need
    df = df[df['AverageAge'].isin(ages)][['AverageAge', 'GMST_50']].set_index('AverageAge').rename(columns={'GMST_50': 'paleoTemp'})

    label = "PhanDA_paleotemps"
    return df, label
    
    


        

