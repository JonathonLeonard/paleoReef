Contained here is the data used in the publication "Coral thermal tolerances constrain extreme Meso-Cenozoic hothouse reconstructions". 

Any questions can be directed to jonathon.leonard@sydney.edu.au

# How to use
- All the stuff you need is in 'plasimGenieWorkflows'. This was initially forked from the paleoreef repository (Salles et al 2025) and for reference those workflows were retained, but I point the user to the original paleoreef directory for that.
- Use the 'preprocessing.ipynb' notebook to get the coral probability maps from the climate data (however pre-computed maps are also included here). This notebook also calculates coral probability from the biomod species distribution model, although we decided against using this method in the main paper (the code is also not well tested so I can't guarantee it will work). 
- The 'Fig-climateScenario.ipynb' notebook can be used to reproduce all the figures produced in the manuscript
- There is also all the python code that contains the classes and functions to do a lot of the nitty gritty stuff

Note that any updates will be to the github (https://github.com/JonathonLeonard/paleoReef) unless there is a major version change.

# Description of data
There is both original and non-original data included here. Below is a description with citations where necessary. If there's no citation, please cite the manuscript attached to this repository.

## biomodEnvRasters
These are all the environmental variable rasters as geoTIFFs so they can be fed into biomod.
## climate_simulations
The PLASIM-GENIE climate simulation outputs.
## cpt_files
Colour palette files created and used for the figures in this manuscript.
## fossils
The database of fossil data for scleractibia corals (from PBDB – https://paleobiodb.org) and the reef data from PARED (Kiessling, Wolfgang, and Cristina Krause. “PaleoReefs Database (PARED).” Version 1.0. Zenodo, February 10, 2022. https://doi.org/10.5281/zenodo.6037852.)
## fosterCO2.txt
Text file that has the median CO2 value by age according to Foster, Gavin L., Dana L. Royer, and Daniel J. Lunt. “Future Climate Forcing Potentially without Precedent in the Last 420 Million Years.” Nature Communications 8, no. 1 (2017): 1. https://doi.org/10.1038/ncomms14845.
## fuzzy_logic_projection
This is all the coral suitability maps based on the transfer function of SST thresholds. It is called fuzzy logic because it began as an alternative to biomod kinda based of a single variable fuzzy logic formula. There is data for both paleogeographies and for each CO2 value, and also the closest co2 value that matches the temperature from the paleotemp curves.
## modelData.xlsx
A bunch of data summaries produced as a supplement to the manuscript attached to this repository.
## paleogeography
The paleoelevation models used in this study:
- 'Scotese' is the PALEOMAP elevations (Scotese, Christopher R, and Nicky Wright. PALEOMAP Paleodigital Elevation Models (PaleoDEMS) for the Phanerozoic. 2018. https://www.earthbyte.org/paleodem-resource-scotese-and-wright-2018/.)
- 'topos_filled_PaleomagRef' is the ZAH22 paleogeography (Zahirovic, Sabin, Ahmed Eleish, Sebastiano Doss, et al. “Subduction and Carbonate Platform Interactions.” Geoscience Data Journal 9, no. 2 (2022): 371–83. https://doi.org/10.1002/gdj3.146.)
 (Cao, Wenchao, Sabin Zahirovic, Nicolas Flament, Simon Williams, Jan Golonka, and R. Dietmar Müller. “Improving Global Paleogeography since the Late Paleozoic Using Paleobiology.” Biogeosciences 14, no. 23 (2017): 5425–39. https://doi.org/10.5194/bg-14-5425-2017.)
## paleotemperature_datasets
The other paleotemperature curve that we compared our coral-constrained curve to:
- Merdith2024 is actually MER25 (Merdith, Andrew S., Thomas M. Gernon, Pierre Maffre, et al. “Phanerozoic Icehouse Climates as the Result of Multiple Solid-Earth Cooling Mechanisms.” Science Advances 11, no. 7 (2025): eadm9798. https://doi.org/10.1126/sciadv.adm9798.)
- paleotemperatures_Grossman_Joachimski_2022 is GRO22 (Grossman, Ethan L., and Michael M. Joachimski. “Ocean Temperatures through the Phanerozoic Reassessed.” Scientific Reports 12, no. 1 (2022): 8938. https://doi.org/10.1038/s41598-022-11493-1.)
- PhanDA_GMSTandCO2_percentiles.csv is the GMST curve from PhanDA (Judd, Emily J., Jessica E. Tierney, Daniel J. Lunt, et al. “A 485-Million-Year History of Earth’s Surface Temperature.” Science 385, no. 6715 (2024): eadk3705. https://doi.org/10.1126/science.adk3705.)
- PhanDA_TropicalSST_30S30N.xlsx is the above but for tropical SST
- scotese2021_paleotemps is SCO21 (Scotese, Christopher R., Haijun Song, Benjamin J. W. Mills, and Douwe G. van der Meer. “Phanerozoic Paleotemperatures: The Earth’s Changing Climate during the Last 540 Million Years.” Earth-Science Reviews 215 (April 2021): 103503. https://doi.org/10.1016/j.earscirev.2021.103503.)
- Song_et_al_2019.xlsx is SONG19 (Song, Haijun, Paul B. Wignall, Huyue Song, Xu Dai, and Daoliang Chu. “Seawater Temperature and Dissolved Oxygen over the Past 500 Million Years.” Journal of Earth Science 30, no. 2 (2019): 236–43. https://doi.org/10.1007/s12583-018-1002-2.)
## photozoan_locations
Is data for the photozoan T factories used to check the apparent reef coral preferred SST ranges.Laugié, Marie, Julien Michel, Alexandre Pohl, Emmanuelle Poli, and Jean Borgomano. “Global Distribution of Modern Shallow-Water Marine Carbonate Factories: A Spatial Model Based on Environmental Parameters.” Scientific Reports 9, no. 1 (2019): 16432. https://doi.org/10.1038/s41598-019-52821-2.
## reconstructions
The tectonic plate reconstructions used here.
- ZAH22 (Zahirovic, Sabin, Ahmed Eleish, Sebastiano Doss, et al. “Subduction and Carbonate Platform Interactions.” Geoscience Data Journal 9, no. 2 (2022): 371–83. https://doi.org/10.1002/gdj3.146.)
- For the PALEOMAP one, use the one included with GPlates (see filepath in notebook)
## sst.mon.ltm.1981-2010.nc
The present day sst map (Huang, B., E. Freeman, J. H. Lawrimore, et al. “Extended Reconstructed Sea Surface Temperature Version 4 (ERSST.v4), Part I. Upgrades and Intercomparisons.” 2014 (December 2014): GC51D-0447. https://ui.adsabs.harvard.edu/abs/2014AGUFMGC51D0447H.)

