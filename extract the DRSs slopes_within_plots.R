# extract the slopes of the diversity-stability relationships (DSRs) at alpha, beta1, and gamma scales
# Load packages
library(nlme)
library(lme4)

# load data
NEON_stab_plots.data <- read.csv('NEON_stab_within_plots.csv',header=T)
# check variables
variable.names(NEON_stab_plots.data)
##############################################################################################################
[1] "no."          "domainID"     "domain_name"  "siteID"       "site_name"    "state"        "MAT_C"       
[8] "MAP_mm"       "Temp_s"       "Prec_s"       "started_year" "end_year"     "duration"     "plotID"      
[15] "latitude"     "longitude"    "elevation_m"  "nlcdClass"    "alpha_div"    "beta_div1"    "gamma_div"   
[22] "spe_sta"      "spe_asy"      "alpha_sta"    "spa_asy1"     "gamma_sta" 
##############################################################################################################

# 1. the relationships of alpha stability to alpha diversity
dd <- NEON_stab_plots.data
ss <- unique(dd$siteID); ss

all_alpha.result <- c()
for (i in 1:length (ss)){
  dd.i <- dd[dd$siteID==ss[i], ]
      each.depend <- dd.i[ , 24] # alpha stability
      each.x <- dd.i[ , 19] # alpha diversity
      lat <- dd.i[ , 15] # latitude
      long <- dd.i[ , 16] # longitude
      each.fit <- gls(each.depend ~ each.x, correlation=corExp(form = ~ lat + long))
      each.result <- c(unlist(summary(each.fit)$tTable[1, 1:2]),
                       unlist(summary(each.fit)$tTable[2, 1:4]))
      d <- c(as.character(ss[i]), each.result)
      all_alpha.result <- rbind(all_alpha.result, d)
}
colnames(all_alpha.result) <- c("siteID","intercept", "se_intercept","slope", "se_slope", "t value", "p")
#write.csv(all_alpha.result, file = "all_alpha.result.csv")


# 2. the relationships of spatial asynchrony to beta diversity across quadrats
dd <- NEON_stab_plots.data
ss <- unique(dd$siteID); ss

all_beta1.result <- c()
for (i in 1:length (ss)){
  dd.i <- dd[dd$siteID==ss[i], ]
  each.depend <- dd.i[ , 25] # spatial asynchrony
  each.x <- dd.i[ , 20] # beta diversity 1
  lat <- dd.i[ , 15] # latitude
  long <- dd.i[ , 16] # longitude
  each.fit <- gls(each.depend ~ each.x, correlation=corExp(form = ~ lat + long))
  each.result <- c(unlist(summary(each.fit)$tTable[1, 1:2]),
                   unlist(summary(each.fit)$tTable[2, 1:4]))
  d <- c(as.character(ss[i]), each.result)
  all_beta1.result <- rbind(all_beta1.result, d)
}
colnames(all_beta1.result) <- c("siteID","intercept", "se_intercept","slope", "se_slope", "t value", "p")
#write.csv(all_beta1.result, file = "all_beta1.result.csv")


# 3. the relationships of gamma stability to gamma diversity 
dd <- NEON_stab_plots.data
ss <- unique(dd$siteID); ss

all_gamma.result <- c()
for (i in 1:length (ss)){
  dd.i <- dd[dd$siteID==ss[i], ]
  each.depend <- dd.i[ , 26] # gamma stability
  each.x <- dd.i[ , 21] # gamma diversity
  lat <- dd.i[ , 15] # latitude
  long <- dd.i[ , 16] # longitude
  each.fit <- gls(each.depend ~ each.x, correlation=corExp(form = ~ lat + long))
  each.result <- c(unlist(summary(each.fit)$tTable[1, 1:2]),
                   unlist(summary(each.fit)$tTable[2, 1:4]))
  d <- c(as.character(ss[i]), each.result)
  all_gamma.result <- rbind(all_gamma.result, d)
}
colnames(all_gamma.result) <- c("siteID","intercept", "se_intercept","slope", "se_slope", "t value", "p")
#write.csv(all_gamma.result, file = "all_gamma.result.csv")
