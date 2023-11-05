library(terra)
library(biomod2)
library(doParallel)
registerDoParallel(cores=8)

res <- 'res2' 
reef <- 'photozoan'
reefName <- paste0(reef,'_',res)

# Read environmental dataset
files <- list.files(path=paste0('env_var/bio_',res,'_0Ma_foster'), pattern='.tif', full.names=TRUE )
nb_env <- length(files)
enviVar <- rast(files[1:nb_env])
crs(enviVar) <- '+proj=longlat +datum=WGS84 +no_defs'
enviVar <- rev(enviVar)
enviVar <- rev(enviVar)

# Read PA file
filename <- file.path(paste0('env_var/factoryReef/',res,'_',reef,'.csv'))
reefFile  <- read.csv(filename)
dfReef <- as.numeric(reefFile[, reef])
reefXY <- reefFile[, c('lon','lat')]

reefBiomodData <- BIOMOD_FormatingData(resp.var = dfReef,
                                     expl.var = enviVar,
                                     resp.xy = reefXY,
                                     resp.name = reefName)

# Default modelling options
reefBiomodOptions <- BIOMOD_ModelingOptions()

# Single models
reefBiomodModelOut <- BIOMOD_Modeling(bm.format = reefBiomodData,
                                    bm.options = reefBiomodOptions,
                                    modeling.id = 'AllModels',
                                    CV.strategy = 'random',
                                    CV.nb.rep = 2,
                                    CV.perc = 0.8,
                                    var.import = 3,
                                    metric.eval = c('TSS','ROC','KAPPA','SR','ACCURACY'),
                                    seed.val = 123,
                                    nb.cpu = 8)
# file.out <- paste0(reefName, "/", reefName, ".AllModels.models.out")
# reefBiomodModelOut <- get(load(file.out))
# reefBiomodModelOut

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