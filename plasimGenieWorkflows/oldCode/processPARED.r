library("palaeoverse")
citation("palaeoverse")

setwd("/Users/jono/OneDrive - The University of Sydney (Staff)/paleoReef/plasimGenieWorkflows")
getwd()

# Load the dataset
data(reefs)
# View the first five rows & columns of reefs
reefs[1:5, 1:5]

# # Call documentation for dataset
# ?reefs
# # You can also use
# help(reefs)

# How many reefs are there in the dataset?
# Each row represents an individual reef.
nrow(reefs)

# How many reefs per interval?
# Let's group by the interval column to test
reef_counts <- group_apply(occdf = reefs, group = "interval", fun = nrow)
head(reef_counts)

# Load the interval key
data("interval_key")
# Assign a common time scale based on an interval key
reefs <- look_up(occdf = reefs,
                 early_interval = "interval",
                 late_interval = "interval",
                 int_key = interval_key)
# Example output
reefs[100:110, 15:19]

write.csv(reefs,'data/PARED_time_simple.csv')

# Now we have numeric ages for our data, we can easily
# remove pre-Phanerozoic data to focus our study
reefs <- subset(reefs, interval_max_ma <= 541)
# Extract Phanerozoic stage-level stages for time bins
bins <- time_bins(interval = "Phanerozoic", rank = "stage")
# Bin data
# bin_time requires "max_ma" and "min_ma" columns
# Rename columns in reefs
colnames(reefs)[which(colnames(reefs) == "interval_max_ma")] <- "max_ma"
colnames(reefs)[which(colnames(reefs) == "interval_min_ma")] <- "min_ma"

# Five methods exist in bin_time for binning occurrence data
# You can see details on each via ?bin_time
# Bin by midpoint age
bin_time(occdf = reefs, bins = bins, method = "mid")
# Bin by overlap majority
# bin_time(occdf = reefs, bins = bins, method = "majority")
# # Bin into every bin within age range
# bin_time(occdf = reefs, bins = bins, method = "all")
# # Bin randomly into bins within age range (equal probability)
# bin_time(occdf = reefs, bins = bins, method = "random", reps = 10)
# # Bin randomly into bins within age range (uniform probability)
# bin_time(occdf = reefs, bins = bins, method = "point", reps = 10)

# Let's go with "all" for this example!
reefs <- bin_time(occdf = reefs, bins = bins, method = "all")

# Count number of occurrences per interval
reefs_time <- group_apply(occdf = reefs, group = "interval_mid_ma", fun = nrow)
# Check output
head(reefs_time)

# install.packages("repr")
library("repr")
options(repr.plot.width=13, repr.plot.height=8)
# Plot data
plot(x = reefs_time$interval_mid_ma,
     y = reefs_time$nrow,
     xlab = "Time (Ma)",
     ylab = "Number of reefs per stage",
     xlim = c(541, 0),
     xaxt = "n",
     type = "b", pch = 20,
#      cex.main=1.25, cex.lab=2, cex.axis=2
    )
# Add axis of the geological time scale
axis_geo(side = 1, intervals = "periods")

# Palaeorotate occurrences
reefs <- palaeorotate(occdf = reefs, age = "bin_midpoint",
                      method = "point", model = "PALEOMAP")

# Check palaeocoordinates
head(reefs[, c("p_lng", "p_lat")])

# Plot data
plot(x = reefs$bin_midpoint,
     y = reefs$p_lat,
     xlab = "Time (Ma)",
     ylab = "Palaeolatitude (\u00B0)",
     xlim = c(541, 0),
     xaxt = "n",
     type = "p", pch = 20)
axis_geo(side = 1, intervals = "periods")

# Let's first assume hemispheric symmetry and convert
# palaeolatitudes to absolute palaeolatitudes
reefs$p_lat <- abs(reefs$p_lat)
# Now we can calculate the most poleward latitude per time bin
# Extract unique interval midpoints
midpoints <- sort(unique(reefs$bin_midpoint))
# Calculate the maximum palaeolatitude for each time bin
reefs_max <- sapply(X = midpoints, FUN = function(x) {
  max(reefs[which(reefs$bin_midpoint == x), ]$p_lat, na.rm = TRUE)
} )
# Plot data
plot(x = midpoints,
     y = reefs_max,
     xlab = "Time (Ma)",
     ylab = "Palaeolatitude (\u00B0)",
     xlim = c(541, 0),
     xaxt = "n",
     type = "b",
     pch = 20)
# Add axis of the geological time scale
axis_geo(side = 1, intervals = "periods")

# We can call multiple models at once with palaeorotate
# First, let's define the models...
models <- c("GOLONKA", "PALEOMAP", "MERDITH2021")
# And now palaeorotate!
reefs <- palaeorotate(occdf = reefs, age = "bin_midpoint",
                      method = "point", model = models)

# Check palaeocoordinates
# When multiple models are called, the name of the model
# is added as a suffix to p_lng and p_lat
head(reefs[, c("p_lng_PALEOMAP", "p_lat_PALEOMAP",
               "p_lng_GOLONKA", "p_lat_GOLONKA",
               "p_lng_MERDITH2021", "p_lat_MERDITH2021")])

write.csv(reefs,'data/PARED_time.csv')