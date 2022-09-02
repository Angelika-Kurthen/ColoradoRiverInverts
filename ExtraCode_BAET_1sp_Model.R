#----------------------------------------\
# Old/Extra Script for BAET_1sp_Model.R
#---------------------------------------



# proportional (linear) relationship between temp and growth/development (maximum parsimony)
#if the water temp is warmer than the average water temp, then development is favored over growth
#P1 = A[1,1] = prob of staying in stage 1
#G1 = A[2,1] = prob of going to stage 2
#P2 = A[2,2] = prob of staying in stage 2
#G1 = A[3,2] = prob of going to stage 3
# if (temps$Temperature[t-1] > mean_temp){
#   ABAET[3,2]= ABAET[2,1] = ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[3,2]) + ABAET[3,2]
#   ABAET[1,1] = ABAET[2,2] = ABAET[2,2] - ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[2,2])
# }
# 
# # if the water temp is cooler than the average water temp, then growth is favored over development
# if (temps$Temperature[t-1] < mean_temp){
#   ABAET[2,2] = ABAET [1,1] = ((mean_temp - temps$Temperature[t-1])/mean_temp * ABAET[2,2]) + ABAET[2,2]
#   ABAET[3,2]= ABAET[2,1] = ABAET[3,2] - ((temps$Temperature[t-1] - mean_temp)/mean_temp * ABAET[3,2])
# }

# 


# using two y = sqrt(x) to model change
#ABAET[3,2] = ABAET[2,1] = 0.656*sqrt(temps$Temperature[t-1]) - 1.79
#ABAET[2,2] = ABAET[1,1] = -0.4051*sqrt(temps$Temperature[t-1]) + 1.62



# regular old logistic equation - issue, can go negative
#if (Total.N[t-1] < K) {       # we can also just imagine 'normal' logistic growth of the K - N/K
#F_BAET <- F_BAET*(1 - Total.N[t-1]/K)
#} else {
# F_BAET <- F_BAET*((1 - (K-1))/K)
#}
