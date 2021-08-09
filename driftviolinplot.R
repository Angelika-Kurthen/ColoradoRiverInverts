
install.packages("devtools")
library(devtools)
install_github(repo = "jmuehlbauer-usgs/R-packages", subdir = "foodbase")
library(foodbase)
library(tidyr)
library(ggplot2)
drift <- readDB(gear = "Drift", type = "Sample", updater = T)


t1 <- readDB(type = 'SpeciesList')

# data from the upper thirds of the river

# NZMS = New Zealand Mudsnail
# GAMM = Gammarus lacustris
# 
# SIML = Simuliidae
# 
species <- c("NZMS","GAMM", "QUAG", "CHIA", "CHIP", "CHIL", "SIMA", "SIMP", "SIML", "HYDE", "BAEL")

for (sp in 1:length(species)){
upper <- drift[drift$RiverMile <= 75, ]
upper <- sampspec(samp = upper)

# data from lower third of the river
lower <- drift[drift$RiverMile >= 150, ]
lower <- sampspec(samp = lower)

# pull species specific data
upper <- sampspec(samp = upper, species = species[sp], stats = T)
# get species data
upper_sp <-  upper$SpecDel
# get rid of NAs
upper_sp <- upper_sp[-which((is.na(upper_sp$CountTotal == "NA"))== T),] 

# make columns rows and rows columns
upper_sp <- gather(upper_sp[ , 3:17], key = "Size")

# remove counts of 0 (since they don't need to be in out distribution)
upper_sp <- upper_sp[which(upper_sp$value != 0), ]

size_bin <- c("B0", "B1", "B2", "B3", "B4", "B5",
              "B6", "B7", "B8","B9", "B10", "B11",
              "B12", "B13")
size_list <- seq(from = 0, to =13)
for(n in 1:length(size_list)){
  upper_sp$Size <- replace(upper_sp$Size, which(upper_sp$Size == size_bin[n]), size_list[n])
  
}

size <- vector()
for (i in 1:length(upper_sp$Size)){
  a <- as.vector(rep(upper_sp$Size[i], times = upper_sp$value[i]))
  size <- c(size, a)
}

location <- rep("RM -15 to 75", length(size))
df <- as.data.frame(cbind(size, location))



# lower third
# pull species specific data
lower <- sampspec(samp = lower, species = "NZMS", stats = T)
# get species data
lower_sp <-  lower$SpecDel
# get rid of NAs
lower_sp <- lower_sp[-which((is.na(lower_sp$CountTotal == "NA"))== T),] 

# make columns rows and rows columns
lower_sp <- gather(lower_sp[ , 3:17], key = "Size")

# remove counts of 0 (since they don't need to be in out distribution)
lower_sp <- lower_sp[which(lower_sp$value != 0), ]

size_bin <- c("B0", "B1", "B2", "B3", "B4", "B5",
              "B6", "B7", "B8","B9", "B10", "B11",
              "B12", "B13")
size_list <- seq(from = 0, to =13)
for(n in 1:length(size_list)){
  lower_sp$Size <- replace(lower_sp$Size, which(lower_sp$Size == size_bin[n]), size_list[n])
  
}

size <- vector()
for (i in 1:length(lower_sp$Size)){
  a <- as.vector(rep(lower_sp$Size[i], times = lower_sp$value[i]))
  size <- c(size, a)
}

location <- rep("RM 150 to 226", length(size))
df.low <- as.data.frame(cbind(size, location))

df <- as.data.frame(rbind(df, df.low))

df$size <- as.numeric(as.character(df$size))

assign(noquote(paste0(species[sp], "_df")), df)

}

# want size proportion - how to standardize by Volume? will not get a whole number?

# now make different plots

# plot of New Zealand Mudsnail
ggplot(data = NZMS_df, aes(x = location, y = size))+
  geom_violin(trim = F)+
  xlab("Location")+
  ylab("Size (mm)")+
  ggtitle("New Zealand Mudsnail Size Distributions")+
  theme_bw()

# plot of Gammarus lacustris
ggplot(data = GAMM_df, aes(x = location, y = size))+
  geom_violin(trim = F)+
  xlab("Location")+
  ylab("Size (mm)")+
  ggtitle("Gammarus lacustris Size Distributions")+
  theme_bw()

# plot of Quagga Mussels
ggplot(data = QUAG_df, aes(x = location, y = size))+
  geom_violin(trim = F)+
  xlab("Location")+
  ylab("Size (mm)")+
  ylim(0, 7)+
  ggtitle("Quagga Mussel Size Distributions")+
  theme_bw()

# combine Chironomidae data and plot

# add stage-based info the size data
stage <- rep("Adult", time = length(CHIA_df$size))
CHIA_df <- cbind(CHIA_df, stage)

stage <- rep("Larvae", time = length(CHIL_df$size))
CHIL_df <- cbind(CHIL_df, stage)

stage <- rep("Pupae", time = length(CHIP_df$size))
CHIP_df <- cbind(CHIP_df, stage)

# combind all the different stage data into one dataframe
Chiro_df <- rbind(CHIA_df, CHIL_df, CHIP_df)

ggplot(data = Chiro_df, aes(x = location, y = size))+
  geom_violin(aes(color = stage), trim= F)+
  xlab("Location")+
  ylab("Size (mm)")+
  ylim(0, 7)+
  ggtitle("Chironomidae Size Distributions")+
  theme_bw()  

# plot Hydropsychidae larvae
ggplot (data = HYDE_df, aes(x = location, y = size))+
  geom_violin(trim = F)+
  xlab("Location")+
  ylab("Size (mm)")+
  ggtitle("Hydrosychidae Larvae Size Distributions")+
  theme_bw()

# plot Baetidae larvae
ggplot(data = BAEL_df, aes(x = location, y= size))+
  geom_violin(trim = F)+
  xlab("Location")+
  ylab("Size (mm)")+
  ggtitle("Baetidae Larvae Size Distributions")+
  theme_bw()

# combine Simuliidae data and plot
# add stage-based info the size data
stage <- rep("Adult", time = length(SIMA_df$size))
SIMA_df <- cbind(SIMA_df, stage)

stage <- rep("Larvae", time = length(SIML_df$size))
SIML_df <- cbind(SIML_df, stage)

stage <- rep("Pupae", time = length(SIMP_df$size))
SIMP_df <- cbind(SIMP_df, stage)

# combind all the different stage data into one dataframe
Sim_df <- rbind(SIMA_df, SIML_df, SIMP_df)

ggplot(data = Sim_df, aes(x = location, y = size))+
  geom_violin(aes(color = stage), trim= F)+
  xlab("Location")+
  ylab("Size (mm)")+
  ylim(0, 7)+
  ggtitle("Simuliidae Size Distributions")+
  theme_bw()  
