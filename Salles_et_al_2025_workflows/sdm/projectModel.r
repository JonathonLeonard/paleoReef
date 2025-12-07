library(terra)
library(biomod2)
library(doParallel)
registerDoParallel(cores=6)

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = 'res2'
}


time <- strtoi(args[1], base=0L)
res <- args[2] 

bres <- 'res2' 
reef <- 'photozoan'
reefName <- paste0(reef,'.',bres)

file.out <- paste0(reefName, "/", reefName, ".AllModels.models.out")
reefBiomodModel <- get(load(file.out))

file2.out <- paste0(reefName, "/", reefName, ".AllModels.ensemble.models.out")
reefBiomodModel_EM <- get(load(file2.out))

# Read environmental dataset
# files <- list.files(path=paste0('env_var_proj/bio_',res,'_',time,'Ma_foster'), pattern='.tif', full.names=TRUE )
files <- list.files(path=paste0('env_var/bio_',res,'_',time,'Ma_foster'), pattern='.tif', full.names=TRUE )
nb_env <- length(files)
evarTime <- rast(files[1:nb_env])
crs(evarTime) <- '+proj=longlat +datum=WGS84 +no_defs'
evarTime <- rev(evarTime)
evarTime <- rev(evarTime)

# Project single models 
reefBiomodProj <- BIOMOD_Projection(bm.mod = reefBiomodModel,
                                  proj.name = paste0(time,'Ma'),
                                  new.env = evarTime,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE)

# Project ensemble models (from single projections)
reefBiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = reefBiomodModel_EM, 
                                             bm.proj = reefBiomodProj,
                                             models.chosen = 'all',
                                             metric.binary = 'all',
                                             metric.filter = 'all')
