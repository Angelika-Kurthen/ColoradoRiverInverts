##################################################################
##Integration of egg survival based on Kennedy et al., 2016
###################################################################

# Under high hydropeaking, egg survival follows eq. 7 in Kennedy et al., 2016

f <- function(x) {exp(-r*2)*x*r}

# where r is where in the river a species lays eggs (BAET = 0, CHIRO = 0.5, SIMU = 0.8)
# c is a constant related to an egg survival curve (eq 2), in this case the c = 2
# f (not included in this calculation) is 0.033
# and lower and upper bounds depend on where in the river a species will lay its eggs
# with 0 being the shore and 1 being the center of the river
# BAET lay eggs on river edge [0,0.2]
# SIMU lay eggs in open water [0.8, 1]
# CHIRO have a middle strategy [0.4, 0.6]

# high hydropeaking egg modifier for BAET
r = 0
integrate(f, lower = 0, upper = 0.2)
H_BAET = 0.003230344*2

# high hydropeaking egg modifier for CHIRO
r = 0.5
integrate(f, lower = 0.4, upper = 0.6)
H_CHIRO = 0.01615172*2

# high hydropeaking egg modifer for SIMU
r = 0.8
integrate(f, lower = 0.8, upper =1 )
H_SIMU = 0.02677842*2
