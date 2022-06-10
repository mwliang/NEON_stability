# Load packages
library(piecewiseSEM)

# load data
NEON_stab.data <- read.csv('NEON_stab.csv',header=T)
# check variables
variable.names(NEON_stab.data)
##############################################################################################################
[1] "no."                       "domain"                    "domain_name"               "siteID"                   
[5] "site_name"                 "state"                     "latitude"                  "longitude"                
[9] "elevation_m"               "MAT_C"                     "MAP_mm"                    "Temp_s"                   
[13] "Prec_s"                    "started_year"              "end_year"                  "duration"                 
[17] "no.plot"                   "spa_dist"                  "alpha_div"                 "beta_div1"                
[21] "gamma_div"                 "beta_div2"                 "tau_div"                   "spe_sta"                  
[25] "spe_asyn"                  "alpha_sta"                 "spa_asyn1"                 "gamma_sta"                
[29] "spa_asyn2"                 "tau_sta"                   "slope_alpha_sta.alpha_div" "slope_spa_asyn1.beta_div1"
[33] "slope_gamma_sta.gamma_div" "se_alpha"                  "se_beta1"                  "se_gamma" 
##############################################################################################################
##############################################################################################################
# Based on our priori structural equation modeling (SEM) (Extended Data Fig 2), 
# we have carried out several SEMs to take account for potential associations among variables.

# First, we have started a full model with all potential pathway effects.

# Second, after attaining a satisfactory model fit, 
# in which the model has low AIC and p > 0.05 without missing significant pathways 
# using Shipley¡¯s test of d-separation (Lefcheck 2016), we then introduced each climatic factors. 

# Lastly, based on the rule of model selection (lower AIC), we have chosen our final SEM.

##############################################################################################################
##############################################################################################################

############################################
# 1. the main SEM in our study (N = 36 at site level)
# the first full SEM, excluding all exogenous variables
SEM.m <- psem(
  lm(spa_asyn1 ~ beta_div1, data = NEON_stab.data), # empirical relationships
  lm(spa_asyn2 ~ beta_div2, data = NEON_stab.data), # empirical relationships
  
  lm(alpha_sta ~ alpha_div, data = NEON_stab.data), # empirical relationships
  lm(gamma_sta ~ alpha_sta + spa_asyn1, data = NEON_stab.data), # mathematical relationships, r2 = 1.00
  lm(tau_sta ~ gamma_sta + spa_asyn2, data = NEON_stab.data), # mathematical relationships, r2 = 1.00
  
  spa_asyn1 %~~% alpha_sta, # a mathematical relationship
  spa_asyn2 %~~% spa_asyn1, # a mathematical relationship
  beta_div1 %~~% beta_div2, # a empirical relationship
  data = NEON_stab.data
)
# To evaluate the model
summary(SEM.m) # Note that all pathway were significant.
fisherC(SEM.m) # Note that SEM.m was unsaturated (Fisher.c > df) and its significant estimations were very well with p > 0.05.
AIC(SEM.m)
dSep(SEM.m) # Note that SEM.m also did not miss any significant pathways.


# second, based on the empirical bivariate relationships of climatic factors to diversity and stability.

# Since the piecewise SEM decomposes the network by combining multiple linear regressions 
# for each response to generate inferences about the entire SEM (Lefcheck 2016), 
# therefore, we checked the bivariate relationships between all climatic variables 
# to ensure that a linear model was appropriate before we introduced these variables into the SEM. 
# Because the introduced variables may mediate the effects strengths of other climatic variables 
# through a significant bivariate correlation. 

# Specifically, precipitation and temperature are highly correlated to their seasonality among 36 NEON sites:
summary(lm(MAP_mm ~ Prec_s, data = NEON_stab.data)) # P = 0.0134
summary(lm(MAT_C ~ Temp_s, data = NEON_stab.data)) # P < 0.00001
summary(lm(MAP_mm ~ MAT_C, data = NEON_stab.data)) # P = 0.016
# Thus, we expected that these seasonality variables still related to MAP/MAT, instead had direct effects on diversity and stability
# We introduced these environmental variables and other abiotic factors in our SEM.m.
SEM_env1.m <- psem(
  lm(alpha_div ~  MAP_mm + MAT_C, data = NEON_stab.data), # empirical relationships
  lm(beta_div1 ~ MAP_mm + MAT_C, data = NEON_stab.data), # empirical relationships
  lm(beta_div2 ~ MAP_mm + MAT_C + no.plot, data = NEON_stab.data), # empirical relationships
  
  lm(spa_asyn1 ~ MAP_mm + MAT_C + beta_div1, data = NEON_stab.data), # empirical relationships
  lm(spa_asyn2 ~ MAP_mm + MAT_C + beta_div2, data = NEON_stab.data), # empirical relationships
  
  lm(alpha_sta ~ MAP_mm + MAT_C + alpha_div, data = NEON_stab.data), # empirical relationships
  lm(gamma_sta ~ alpha_sta + spa_asyn1, data = NEON_stab.data), # mathematical relationships, r2 = 1.00
  lm(tau_sta ~ gamma_sta + spa_asyn2, data = NEON_stab.data), # mathematical relationships, r2 = 1.00
  
  spa_asyn1 %~~% alpha_sta, # a mathematical relationship
  spa_asyn2 %~~% spa_asyn1, # a mathematical relationship
  beta_div1 %~~% beta_div2, # a empirical relationship
  MAP_mm %~~% MAT_C, # a empirical relationship from above significant bivariate correlation 
  MAP_mm %~~% Prec_s, # a empirical relationship from above significant bivariate correlation 
  MAT_C %~~% Temp_s, # a empirical relationship from above significant bivariate correlation 
  data = NEON_stab.data
)

# To evaluate the model
summary(SEM_env1.m) # Note that there were non-significant pathways.
fisherC(SEM_env1.m) # Note that SEM_env1.m was still unsaturated (Fisher.c > df) and its significant estimations were very well with p > 0.05.
AIC(SEM_env1.m)
dSep(SEM_env1.m) # Note that SEM_env1.m did have a missing significant pathway.



# And then, we fitted the SEM_env2.m without these non-significant pathways.
SEM_env2.m <- psem(
  #lm(alpha_div ~  MAP_mm + MAT_C, data = NEON_stab.data), # empirical relationships
  lm(beta_div1 ~ MAP_mm + MAT_C, data = NEON_stab.data), # empirical relationships
  lm(beta_div2 ~ MAT_C + no.plot, data = NEON_stab.data), # empirical relationships
  
  lm(spa_asyn1 ~ beta_div1, data = NEON_stab.data), # empirical relationships
  lm(spa_asyn2 ~ MAP_mm + beta_div2, data = NEON_stab.data), # empirical relationships
  
  lm(alpha_sta ~ MAT_C + alpha_div, data = NEON_stab.data), # empirical relationships
  lm(gamma_sta ~ alpha_sta + spa_asyn1, data = NEON_stab.data), # mathematical relationships, r2 = 1.00
  lm(tau_sta ~ gamma_sta + spa_asyn2, data = NEON_stab.data), # mathematical relationships, r2 = 1.00
  
  spa_asyn1 %~~% alpha_sta, # a mathematical relationship
  spa_asyn2 %~~% spa_asyn1, # a mathematical relationship
  beta_div1 %~~% beta_div2, # a empirical relationship
  MAP_mm %~~% MAT_C, # a empirical relationship from above significant bivariate correlation 
  MAP_mm %~~% Prec_s, # a empirical relationship from above significant bivariate correlation 
  MAT_C %~~% Temp_s, # a empirical relationship from above significant bivariate correlation 
  data = NEON_stab.data
)

# To evaluate the model
summary(SEM_env2.m) # Note that all pathway were significant.
fisherC(SEM_env2.m) # Note that SEM_env2.m was still unsaturated (Fisher.c > df) and its significant estimations were very well with p > 0.05.
AIC(SEM_env2.m)
dSep(SEM_env2.m) # Note that SEM_env2.m did not have any missing significant pathway.



# Since both Prec_s and Temp_s had no significant effects on diversity and stability, expect for correlations with MAP and MAT, respectively, 
# we removed both Prec_s and Temp_s from SEM_env2.m to simplify our model.
SEM_env3.m <- psem(
  #lm(alpha_div ~  MAP_mm + MAT_C, data = NEON_stab.data), # empirical relationships
  lm(beta_div1 ~ MAP_mm + MAT_C, data = NEON_stab.data), # empirical relationships
  lm(beta_div2 ~ MAT_C + no.plot, data = NEON_stab.data), # empirical relationships
  
  lm(spa_asyn1 ~ beta_div1, data = NEON_stab.data), # empirical relationships
  lm(spa_asyn2 ~ MAP_mm + beta_div2, data = NEON_stab.data), # empirical relationships
  
  lm(alpha_sta ~ MAT_C + alpha_div, data = NEON_stab.data), # empirical relationships
  lm(gamma_sta ~ alpha_sta + spa_asyn1, data = NEON_stab.data), # mathematical relationships, r2 = 1.00
  lm(tau_sta ~ gamma_sta + spa_asyn2, data = NEON_stab.data), # mathematical relationships, r2 = 1.00
  
  spa_asyn1 %~~% alpha_sta, # a mathematical relationship
  spa_asyn2 %~~% spa_asyn1, # a mathematical relationship
  beta_div1 %~~% beta_div2, # a empirical relationship
  MAP_mm %~~% MAT_C, # a empirical relationship from above significant bivariate correlation 
  #MAP_mm %~~% Prec_s, # a empirical relationship from above significant bivariate correlation 
  #MAT_C %~~% Temp_s, # a empirical relationship from above significant bivariate correlation 
  data = NEON_stab.data
)

# To evaluate the model
summary(SEM_env3.m) # Note that all pathway were significant.
fisherC(SEM_env3.m) # Note that SEM_env3.m was still unsaturated (Fisher.c > df) and its significant estimations were very well with p > 0.05.
AIC(SEM_env3.m) 
dSep(SEM_env3.m) # Note that SEM_env3.m did not have any missing significant pathway.

# Note that all estimations of SEM_env3.m (Std.Estimate, Fisher.c, p, and AIC) were consistent with that of SEM_env2.m.
# We can chose SEM_env3.m as our finial SEM!

##############################################################################################################
##############################################################################################################
# 2. the sub_SEM in our study (N = 945 at plot level)
# This sub-dataset was hierarchical at plot level within 36 NEON sites, 
# within different vegetation types (e.g., the National Land Cover Database Vegetation Type Name, nlcdClass).
# Thus, we have chosen the mixed-effects models, with "siteID/nlcdClass" as a random factor, 
# including "correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass)" as a spatial correlation structure.
 
# Load packages
library(nlme)

# load data
NEON_stab_plots.data <- read.csv('NEON_stab_within_plots.csv',header=T)
# check variables
variable.names(NEON_stab_plots.data)
##############################################################################################################
[1] "no."          "domainID"     "domain_name"  "siteID"       "site_name"    "state"        "MAT_C"       
[8] "MAP_mm"       "Temp_s"       "Prec_s"       "started_year" "end_year"     "duration"     "plotID"      
[15] "latitude"     "longitude"    "elevation_m"  "nlcdClass"    "no.plots"     "alpha_div"    "beta_div1"   
[22] "gamma_div"    "spe_sta"      "spe_asy"      "alpha_sta"    "spa_asy1"     "gamma_sta" 
##############################################################################################################

# the first full sub_SEM, excluding all exogenous variables
sub_SEM.m <- psem(
  
  lme(alpha_sta ~ alpha_div, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # a empirical relationship
  lme(spa_asy1 ~ beta_div1, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # a empirical relationship
  
  lm(gamma_sta ~ alpha_sta + spa_asy1, data = NEON_stab_plots.data), # mathematical relationships, r2 = 1.00
  
  alpha_sta %~~% spa_asy1, # a mathematical relationship
  alpha_div %~~% beta_div1, # a empirical relationship
  data = NEON_stab_plots.data
)
# To evaluate the model
summary(sub_SEM.m) # Note that all pathway were significant.
fisherC(sub_SEM.m) # Note that sub_SEM.m was unsaturated (Fisher.c > df) and its significant estimations were very well with p > 0.05.
AIC(sub_SEM.m)
dSep(sub_SEM.m) # Note that sub_SEM.m had a missing significant pathway.



# Based on the model fit information, we added the significant pathway
sub_SEM_1.m <- psem(
  
  lme(alpha_sta ~ alpha_div + beta_div1, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # a empirical relationship
  lme(spa_asy1 ~ beta_div1, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # a empirical relationship
  
  lm(gamma_sta ~ alpha_sta + spa_asy1, data = NEON_stab_plots.data), # mathematical relationships, r2 = 1.00
  
  alpha_sta %~~% spa_asy1, # a mathematical relationship
  alpha_div %~~% beta_div1, # a empirical relationship
  data = NEON_stab_plots.data
)
# To evaluate the model
summary(sub_SEM_1.m) # Note that all pathway were significant.
fisherC(sub_SEM_1.m) # Note that sub_SEM_1.m was saturated (Fisher.c < df) and its significant estimations were very well with p > 0.05.
AIC(sub_SEM_1.m)
dSep(sub_SEM_1.m) # Note that sub_SEM.m did not have any missing significant pathway.



# Following up the above process and results, we introduced MAP and MAT into sub_SEM_1.m.
sub_SEM_2.m <- psem(
  
  lme(alpha_div ~ MAP_mm + MAT_C, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # empirical relationships
  lme(beta_div1 ~ MAP_mm + MAT_C, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # empirical relationships
  
  lme(alpha_sta ~ alpha_div + beta_div1, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # a empirical relationship
  lme(spa_asy1 ~ beta_div1, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # a empirical relationship
  
  lm(gamma_sta ~ alpha_sta + spa_asy1, data = NEON_stab_plots.data), # mathematical relationships, r2 = 1.00
  
  alpha_sta %~~% spa_asy1, # a mathematical relationship
  alpha_div %~~% beta_div1, # a empirical relationship
  MAP_mm %~~%  MAT_C,  # a mathematical relationship
  data = NEON_stab_plots.data
)
# To evaluate the model
summary(sub_SEM_2.m) # Note that there were two non-significant pathways.
fisherC(sub_SEM_2.m) # Note that sub_SEM_2.m was saturated (Fisher.c < df) and its significant estimations were very well with p > 0.05.
AIC(sub_SEM_2.m)
dSep(sub_SEM_2.m) # Note that sub_SEM_2.m had a missing significant pathway.



# And then, we removed the non-significant pathway and added the significant pathway based one the model fit information.
sub_SEM_3.m <- psem(
  
  #lme(alpha_div ~ MAP_mm + MAT_C, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # empirical relationships
  lme(beta_div1 ~ MAP_mm + MAT_C, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # empirical relationships
  
  lme(alpha_sta ~ alpha_div + beta_div1 + MAT_C, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # a empirical relationship
  lme(spa_asy1 ~ beta_div1, random = ~1|siteID/nlcdClass, correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots.data), # a empirical relationship
  
  lm(gamma_sta ~ alpha_sta + spa_asy1, data = NEON_stab_plots.data), # mathematical relationships, r2 = 1.00
  
  alpha_sta %~~% spa_asy1, # a mathematical relationship
  alpha_div %~~% beta_div1, # a empirical relationship
  MAP_mm %~~%  MAT_C,  # a mathematical relationship
  data = NEON_stab_plots.data
)
# To evaluate the model
summary(sub_SEM_3.m) # Note that all pathway were significant.
fisherC(sub_SEM_3.m) # Note that sub_SEM_3.m was saturated (Fisher.c < df) and its significant estimations were very well with p > 0.05.
AIC(sub_SEM_3.m)
dSep(sub_SEM_3.m) # Note that sub_SEM.m did not have any missing significant pathway.
# We can chose SEM_env3.m as our finial sub_SEM!

