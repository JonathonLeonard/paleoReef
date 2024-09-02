library(terra)
library(biomod2)
library(doParallel)
registerDoParallel(cores=6)

setwd("/Users/jono/OneDrive - The University of Sydney (Staff)/paleoReef/presentDayModelSetup")
getwd()

list.files(path='data')

photo <- read.csv("data/photozoan_filtered.csv")
head(photo)
summary(photo)

# Select the name of the studied species
respName <- 'photozoan'

photozoan <- as.numeric(photo[, respName])

# Get corresponding XY coordinates
photozoanXY <- photo[, c('lon', 'lat')]

# Read environmental dataset
files <- list.files(path=paste0('data/rasters/0Ma_envRaster'), pattern='.tif', full.names=TRUE )
nb_env <- length(files)
enviVar <- rast(files[1:nb_env])
crs(enviVar) <- '+proj=longlat +datum=WGS84 +no_defs'
enviVar <- rev(enviVar)
enviVar <- rev(enviVar)

# bioclim_ZA_sub <- 
#   raster::stack(
#     c(
#       bio_5  = "data/raster/SST.tif",
#       bio_7  = 'data/rasters/SSS.tif',
#       bio_11 = 'data/rasters/netsolar.tif'
#     #   bio_19 = '../data/worldclim_ZA/netsolar.tif'
#     )
#   )

reefBiomodData <- BIOMOD_FormatingData(resp.var = photozoan,
                                     expl.var = enviVar,
                                     resp.xy = photozoanXY,
                                     resp.name = respName,
                                    PA.nb.rep = 2,
                                     PA.nb.absences = 500,
                                    PA.strategy = 'random',
                                    filter.raster = TRUE)

# ## format the data ----
# ProLau_data <- 
#   BIOMOD_FormatingData(
#     resp.var = myResp,
#     resp.xy = myRespXY,
#     expl.var = bioclim_ZA_sub,
#     resp.name = myRespName,
#     PA.nb.rep = 2,
#     PA.nb.absences = 500,
#     PA.strategy = 'random'
#   )

## formatted object summary
reefBiomodData

# chosenModels <- c('GAM', 'GBM', 'GLM')
chosenModels <- c('ANN','FDA','GAM', 'GBM', 'GLM','MARS','RF','SRE','XGBOOST')

# Default modelling options
# reefBiomodOptions <- bm_ModelingOptions(data.type = 'binary',
#                                         strategy = 'default',
#                                         models = chosenModels)

print('test -3')
# Single models
reefBiomodModelOut <- BIOMOD_Modeling(bm.format = reefBiomodData,
                                    # OPT.user = reefBiomodOptions,
                                    models = chosenModels,
                                    modeling.id = 'AllModels',
                                    # OPT.strategy = 'bigboss',
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    var.import = 3,
                                    metric.eval = c('TSS','ROC','KAPPA','SR','ACCURACY'),
                                    seed.val = 123,
                                    nb.cpu = 4)


# Ensemble models
reefBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = reefBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      em.algo = c('EMmean', 'EMcv', 'EMci', 'EMmedian', 'EMca', 'EMwmean'),
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.8),
                                      metric.eval = c('TSS', 'ROC','KAPPA','SR','ACCURACY'),
                                      var.import = 3,
                                      EMci.alpha = 0.05,
                                      EMwmean.decay = 'proportional',
                                      prob.mean = TRUE, 
                                      prob.cv = TRUE, 
                                      prob.ci = TRUE, 
                                      prob.ci.alpha = 0.05, 
                                      prob.median = FALSE, 
                                      committee.averaging = TRUE,
                                      prob.mean.weight = TRUE,
                                      nb.cpu = 8)

# # file.out <- paste0('photozoan/models/AllModels/.AllModels.models.out')
# # reefBiomodModelOut <- get(load(file.out))
# # reefBiomodModelOut

# # Get evaluation scores & variables importance
# get_evaluations(reefBiomodModelOut)
# get_variables_importance(reefBiomodModelOut)

# # Represent evaluation scores & variables importance
# bm_PlotEvalMean(bm.out = reefBiomodModelOut)
# bm_PlotEvalBoxplot(bm.out = reefBiomodModelOut, group.by = c('algo', 'algo'))
# bm_PlotEvalBoxplot(bm.out = reefBiomodModelOut, group.by = c('algo', 'run'))
# bm_PlotVarImpBoxplot(bm.out = reefBiomodModelOut, group.by = c('expl.var', 'algo', 'algo'))
# bm_PlotVarImpBoxplot(bm.out = reefBiomodModelOut, group.by = c('expl.var', 'algo', 'run'))
# bm_PlotVarImpBoxplot(bm.out = reefBiomodModelOut, group.by = c('algo', 'expl.var', 'run'))
# print('test -2')
# # Represent response curves
# bm_PlotResponseCurves(bm.out = reefBiomodModelOut, 
#                       models.chosen = get_built_models(reefBiomodModelOut)[c(1:3, 12:14)],
#                       fixed.var = 'median')
# bm_PlotResponseCurves(bm.out = reefBiomodModelOut, 
#                       models.chosen = get_built_models(reefBiomodModelOut)[c(1:3, 12:14)],
#                       fixed.var = 'min')
# bm_PlotResponseCurves(bm.out = reefBiomodModelOut, 
#                       models.chosen = get_built_models(reefBiomodModelOut)[3],
#                       fixed.var = 'median',
#                       do.bivariate = TRUE)
# print('test -1')
# # Project at 120 Ma
# # Read environmental dataset
# files <- list.files(path=paste0('data/rasters/120Ma_envRaster'), pattern='.tif', full.names=TRUE )
# nb_env <- length(files)
# print(nb_env)
# enviVar <- rast(files[1:nb_env])
# crs(enviVar) <- '+proj=longlat +datum=WGS84 +no_defs'
# enviVar <- rev(enviVar)
# enviVar <- rev(enviVar)


# # Project single models
# reefBiomodProj <- BIOMOD_Projection(bm.mod = reefBiomodModelOut,
#                                   proj.name = '120Ma',
#                                   new.env = enviVar,
#                                   models.chosen = 'all',
#                                   metric.binary = 'all',
#                                   metric.filter = 'all',
#                                   build.clamping.mask = TRUE)

# reefBiomodProj
# plot(reefBiomodProj)
