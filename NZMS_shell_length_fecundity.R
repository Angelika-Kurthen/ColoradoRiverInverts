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