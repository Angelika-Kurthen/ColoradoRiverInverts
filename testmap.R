#===============
# LOAD PACKAGES
#===============
library(tidyverse)
library(maptools)
#===============
# GET RIVER DATA
#===============
#==========
# LOAD DATA
#==========
#DEFINE URL
# - this is the location of the file
url.river_data <- url("http://sharpsightlabs.com/wp-content/datasets/usa_rivers.RData")
# LOAD DATA
# - this will retrieve the data from the URL
load(url.river_data)
# INSPECT
summary(lines.rivers)
lines.rivers@data %>% glimpse()
levels(lines.rivers$FEATURE)
table(lines.rivers$FEATURE)
#==============================================
# REMOVE MISC FEATURES
# - there are some features in the data that we
#   want to remove
#==============================================
lines.rivers <- subset(lines.rivers, !(FEATURE %in% c("Shoreline"
                                                      ,"Shoreline Intermittent"
                                                      ,"Null"
                                                      ,"Closure Line"
                                                      ,"Apparent Limit"
)))
# RE-INSPECT
table(lines.rivers$FEATURE)
#==============
# REMOVE STATES
#==============
#-------------------------------
# IDENTIFY STATES
# - we need to find out
#   which states are in the data
#-------------------------------
table(lines.rivers$STATE)
#---------------------------------------------------------
# REMOVE STATES
# - remove Alaska, Hawaii, Puerto Rico, and Virgin Islands
# - these are hard to plot in a confined window, so 
#   we'll remove them for convenience
#---------------------------------------------------------
lines.rivers <- subset(lines.rivers, !(STATE %in% c('AK','HI','PR','VI')))
# RE-INSPECT
table(lines.rivers$STATE)
#============================================
# FORTIFY
# - fortify will convert the 
#   'SpatialLinesDataFrame' to a proper
#    data frame that we can use with ggplot2
#============================================
df.usa_rivers <- fortify(lines.rivers)
#============
# GET USA MAP
#============
map.usa_country <- map_data("usa")
map.usa_states <- map_data("state")
#=======
# PLOT
#=======
ggplot() +
  geom_polygon(data = map.usa_country, aes(x = long, y = lat, group = group), fill = "#484848") +
  geom_path(data = df.usa_rivers, aes(x = long, y = lat, group = group), color = "#8ca7c0", size = .08) +
  coord_map(projection = "albers", lat0 = 30, lat1 = 40, xlim = c(-121,-73), ylim = c(25,51)) +
  labs(title = "Rivers and waterways of the United States") +
  annotate("text", label = "sharpsightlabs.com", color = "#A1A1A1"
           , x = -89, y = 26.5, size = 5) +
  theme(panel.background = element_rect(fill = "#292929")
        ,plot.background = element_rect(fill = "#292929")
        ,panel.grid = element_blank()
        ,axis.title = element_blank()
        ,axis.text = element_blank()
        ,axis.ticks = element_blank()
        ,text = element_text(color = "#A1A1A1")
        ,plot.title = element_text(size = 34)
  ) 
