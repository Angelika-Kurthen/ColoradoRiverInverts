
install.packages("devtools")
library(devtools)
install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)
drift.data.total <- readDB(gear = "Drift", type = "Sample", updater = TRUE)

# we will say that Lees Ferry = -2 km to 2 km (4k reach)

drift.LF <- drift.data.total[which(drift.data.total$RiverMile >= -2 & drift.data.total$RiverMile <= 2),]

Baet.samp.LF <- sampspec(samp = drift.LF, species = "BAET", stats = T)

