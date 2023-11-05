# Abiotic controls on warm-water carbonates through geological times

Series of workflows used to evaluating global carbonate distribution and accumulation over the past 265 Myr accounting for changes in paleogeography, plate tectonics, and paleo-climatic conditions. 

This is done using a specific carbonate factory, the _photozoan-T factory_, that encompasses both tropical carbonate shelves and detached rimmed platforms and is composed of most of the highest carbonate producers including corals, stromatoporoids, green algae, rudists, and photosymbiotic foraminifers ([more info about the factory](https://www.nature.com/articles/s41598-019-52821-2)).

<img width="1096" alt="SDM_model" src="https://github.com/Geodels/paleoReef/assets/7201912/58ad164e-e4ae-4933-a531-5ce4a4b8cc11">

We rely on a _species distribution modelling_ (SDM) approach to estimate photozoan-T carbonate factory paleo-distribution based on the package [BIOMOD](https://biomodhub.github.io/biomod2/). 

The workflows, presented in the repository in the form of Jupyter Notebooks, describe the different pre- and post-processing scripts that were used to: predict carbonate suitable habitats paleo-distribution, spatio-temporal maps of carbonate accumulations and their associated volumes. 

The approach relies on _(1)_ paleo-environmental ocean variables from the HadCM3 ([BRIDGE](https://www.paleo.bristol.ac.uk)) coupled atmosphere-ocean-vegetation Hadley Centre climate models, _(2)_ the observed modern photozoan-T factory defined after [Michel et al.](https://www.researchgate.net/profile/Julien-Michel-5/publication/333885781_Marine_carbonate_factories_a_global_model_of_carbonate_platform_distribution/links/5d3b098e299bf1995b4cd0ad/Marine-carbonate-factories-a-global-model-of-carbonate-platform-distribution.pdf), _(3)_ the updated continental margins and paleocoastlines dataset available on [Zenodo](https://doi.org/10.5281/zenodo.3903163) based on the paleogeography maps from the [PALEOMAP](https://zenodo.org/records/5460860) project, and _(4)_ sediment flux derived from global scale landscape evolution model for the Phanerozoic [HydroShare](www.hydroshare.org/resource/0106c156507c4861b4cfd404022f9580). In addition, sensitivity tests are making use of the fossil records from the Paleo Reefs PARED [database](https://www.paleo-reefs.pal.uni-erlangen.de). 

