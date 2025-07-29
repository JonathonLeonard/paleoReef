# Predicting paleo-reef distributions

Adaptation of the [workflows developed by Tristan Salles](https://github.com/Geodels/paleoReef) that predict reef distribution using the HADCM3 set of climate simulations. This original workflow is found in the folder *HADCM3Workflows*.

## To run with plasim-genie climate output
This is all done in the *plasimGenieWorkflows* folder. 

1. Get present day carbonate locations in `getCarbonateLocations.ipynb`
2. Get the preferred arrangement of CO<sup>2</sup> in `GMST_CO2.ipynb`
3. Get the environmental rasters that biomod reads in the `getEnvRasters.py` script (not the jupyter notebook). At the moment biomod will read every raster in the folder, so this needs to be re-run everytime you change the set of rasters used by biomod. This is just done through commenting out different raster save lines of code. 
4. Use `processpared.ipynb` to get a csv of all the PARED coral instances that we want to use along with their reconstructed paleo-location. 
5. Use `runModel.ipynb` to run the biomod model, but also to clip the data to shelves, and plot the data. A lot of work is done in this notebook. 

------
*Note*: runModel.ipynb calls  the two R scripts `getEnsembleModel.r` and `getEnsembleProjection.r`