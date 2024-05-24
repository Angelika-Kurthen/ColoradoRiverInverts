##----------------------------------------
## Shell Length and Fecundity relationship
##----------------------------------------

# read in excel supplementary material from McKenzie et al 2017
library(readxl)
data <- read_excel("C:/Users/jelly/Downloads/cjz-2012-0183suppl.xls", col_names = FALSE)

# isolate length and embryo data
length_fec_data <- as.data.frame(data[-c(1:2), c(9,15)])
colnames(length_fec_data) <- c("Length (mm)", "Total Embryos")
length_fec_data$`Length (mm)` <- as.numeric(length_fec_data$`Length (mm)`)
length_fec_data$`Total Embryos` <- as.numeric(length_fec_data$`Total Embryos`) 


lm <- lm(length_fec_data$`Total Embryos` ~ length_fec_data$`Length (mm)`)
summary(lm) # Fecundity = -50.2907 + (16.5408 * x)

# predicted fecundity at 3.2 = 2.63986
-50.2907 + (16.5408 * 3.2)
# predicted fecundity at 3.9538776 = 15.1096
-50.2907 + (16.5408 * 3.9538776)
# predicted fecundity at 5.5 = 40.6837
-50.2907 + (16.5408 * 5.5)
# 


# growth - according to _ growth is dependent on shell length, as is fecundity 
# we want 1 class of non-reproductive subadults (smaller than 3.2 mm) and two classes of reproductive adults (3.2 mm and larger)
# growth follows the equation growth per day = -0.006*length[t]+0.029 (Cross et al., 2010)
# this can be re-written as length[t+1] = 0.994*length[t]+0.029
shell.growth(0.994, 0.029, 0.5)
# on day 164, length ~ 3.2 (maturity) - so on week 24, maturity reached [stage 1 = 0.5mm - 3.2 mm] in about 24 weeks or 12 timesteps
# on day 266, length ~ 3.9538776 - [stage 2 = 3.2 - 3.9538776] in about 14 weeks or 7 timesteps
# that means [stage 3 = 3.9538776+] for about 14 weeks or 7 timesteps
# growth is related to temperature