##################################################################
##Integration of egg survival based on Kennedy et al., 2016
###################################################################

# Under high hydropeaking, egg survival follows eq. 7 in Kennedy et al., 2016

hydrofunction <- function(x) {exp(-x*c)}

# where r is where in the river a species lays eggs (BAET = 0, CHIRO = 0.5, SIMU = 0.8)
# c is a constant related to an egg survival curve (eq 2), in this case the c = 2
# f (not included in this calculation) is 0.033
# and lower and upper bounds depend on where in the river a species will lay its eggs
# with 0 being the shore and 1 being the center of the river
# BAET lay eggs on river edge [0,0.2]
# SIMU lay eggs in open water [0.8, 1]
# CHIRO have a middle strategy [0.4, 0.6]

# high hydropeaking egg modifier for BAET
h = 0.56
# proportion of eggs lost to desiccation
int_baet <- integrate(f, lower = 0, upper = 0.2)
H_BAET = (1 - int_baet$value*2) *(1-h)

# high hydropeaking egg modifier for CHIRO (and HYOS?)
h = 0.56
int_chiro <- integrate(f, lower = 0.4, upper = 0.6)
H_CHIRO = (1 - int_chiro$value*2) *(1-h)

# high hydropeaking egg modifer for SIMU
h = 0.56
int_simu <- integrate(f, lower = 0.8, upper =1 )
H_SIMU = (1- int_simu$value*2) * (1-h)


hydropeaking.mortality <- function(lower, upper, h){
  int <- integrate(hydrofunction, lower = lower, upper = upper)
  H <- (1 - int$value*2) * (1-h)
  return(H)
}
