# Load library
library(piecewiseSEM)

# load data
NEON_stab.data <- read.csv('NEON_stab.csv',header=T)
# check variables
variable.names(NEON_stab.data)
##############################################################################################################
[1] "no."                       "domain"                    "domain_name"              
[4] "siteID"                    "site_name"                 "state"                    
[7] "latitude"                  "longitude"                 "elevation_m"              
[10] "MAT_C"                     "MAP_mm"                    "Temp_s"                   
[13] "Prec_s"                    "started_year"              "end_year"                 
[16] "duration"                  "no.plot"                   "alpha_div"                
[19] "beta1"                     "gamma_div"                 "beta2"                    
[22] "tau_div"                   "spe_sta"                   "spe_asyn"                 
[25] "alpha_sta"                 "spa_asyn1"                 "gamma_sta"                
[28] "spa_asyn2"                 "tau_sta"                   "slope_alpha_sta.alpha_div"
[31] "slope_beta_beta"           "slope_gamma_gamma"  
##############################################################################################################
##############################################################################################################
SEM.m <- psem(
  #lm(alpha_div ~ MAP_mm + MAT_C, data = NEON_stab.data),
  lm(beta1 ~ MAP_mm + MAT_C, data = NEON_stab.data),
  lm(beta2 ~ no.plot  + MAT_C, data = NEON_stab.data),
  
  lm(spa_asyn1 ~ beta1, data = NEON_stab.data),
  lm(spa_asyn2 ~ beta2 + MAP_mm, data = NEON_stab.data),
  
  lm(alpha_sta ~ alpha_div + MAT_C, data = NEON_stab.data),
  lm(gamma_sta ~ alpha_sta + spa_asyn1, data = NEON_stab.data),
  lm(tau_sta ~ gamma_sta + spa_asyn2, data = NEON_stab.data),
  
  spa_asyn1 %~~% alpha_sta,
  spa_asyn2 %~~% spa_asyn1,
  beta1 %~~% beta2,
  MAP_mm %~~% MAT_C,
  data = NEON_stab.data
)
fisherC(SEM.m)
summary(SEM.m)
AIC(SEM.m)
dSep(SEM.m)

