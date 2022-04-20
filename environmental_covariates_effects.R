# Load packages
library(MuMIn)
library(ggplot2)

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

########################################
my_scale <- function(x){
  x1 <- (x-mean(x))/sd(x)
  return(x1)
}
########################################

############################################
# 1 for >= 4-year observations (N=36)
NEON_stab_4yr.data <- NEON_stab.data
data_4yr <- NEON_stab_4yr.data

# 1.1 mean annual precipitation (MAP) (Figs. 4D-4F, Extended Data Fig. 6, Supplementary Tables 8-11)
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(MAP_mm), data = NEON_stab_4yr.data))
summary(lm(tau_div ~ my_scale(MAP_mm), data = NEON_stab_4yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spa_asyn1 ~ my_scale(MAP_mm), data = NEON_stab_4yr.data))
summary(lm(spa_asyn2 ~ my_scale(MAP_mm), data = NEON_stab_4yr.data))
summary(lm(tau_sta ~ my_scale(MAP_mm), data = NEON_stab_4yr.data))

all_MAP_4yr.result <- c()
for (i in 19:33){
  each.depend <- data_4yr[ , i]
  each.x <- data_4yr[ , 11] # MAP
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_4yr)[i], each.result, r2)
  all_MAP_4yr.result <- rbind(all_MAP_4yr.result, d)
}
colnames(all_MAP_4yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_MAP_4yr.result, file = "all_MAP_effects_4yr.csv")

# 1.2 mean annual temperature (MAT) (Figs. 4D-4F, Extended Data Fig. 6, Supplementary Tables 8-11)
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(MAT_C)*my_scale(MAP_mm), data = NEON_stab_4yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spe_sta ~ my_scale(MAT_C), data = NEON_stab_4yr.data))
summary(lm(alpha_sta ~ my_scale(MAT_C), data = NEON_stab_4yr.data))

all_MAT_4yr.result <- c()
for (i in 19:33){
  each.depend <- data_4yr[ , i]
  each.x <- data_4yr[ , 10] # MAT
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_4yr)[i], each.result, r2)
  all_MAT_4yr.result <- rbind(all_MAT_4yr.result, d)
}
colnames(all_MAT_4yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_MAT_4yr.result, file = "all_MAT_effects_4yr.csv")

# 1.3 precipitation seasonality (Prec_s) (Figs. 4D-4F, Extended Data Fig. 6, Supplementary Tables 8-11)
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(Prec_s), data = NEON_stab_4yr.data))
summary(lm(beta_div2 ~ my_scale(Prec_s), data = NEON_stab_4yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spa_asyn1 ~ my_scale(Prec_s), data = NEON_stab_4yr.data))
summary(lm(spa_asyn2 ~ my_scale(Prec_s), data = NEON_stab_4yr.data))
summary(lm(tau_sta ~ my_scale(Prec_s), data = NEON_stab_4yr.data))
# the significant (P < 0.05) effects on the slopes of diversity-stability relationships
summary(lm(slope_alpha_sta.alpha_div ~ my_scale(Prec_s), data = NEON_stab_4yr.data))

all_Prec_s_4yr.result <- c()
for (i in 19:33){
  each.depend <- data_4yr[ , i]
  each.x <- data_4yr[ , 13] #Prec_s
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_4yr)[i], each.result, r2)
  all_Prec_s_4yr.result <- rbind(all_Prec_s_4yr.result, d)
}
colnames(all_Prec_s_4yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_Prec_s_4yr.result, file = "all_Prec_s_effects_4yr.csv")

# 1.4 temperature seasonality (Temp_s) (Figs. 4D-4F, Extended Data Fig. 6, Supplementary Tables 8-11)
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(Temp_s), data = NEON_stab_4yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spe_sta ~ my_scale(Temp_s), data = NEON_stab_4yr.data))

all_Temp_s_4yr.result <- c()
for (i in 19:33){
  each.depend <- data_4yr[ , i]
  each.x <- data_4yr[ , 12] # Temp_s
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_4yr)[i], each.result, r2)
  all_Temp_s_4yr.result <- rbind(all_Temp_s_4yr.result, d)
}
colnames(all_Temp_s.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_Temp_s_4yr.result, file = "all_Temp_s_effects_4yr.csv")

# 1.5 the number of plots (no.plot) (Supplementary Table 12)
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div2 ~ my_scale(no.plot), data = NEON_stab_4yr.data))
summary(lm(tau_div ~ my_scale(no.plot), data = NEON_stab_4yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spe_sta ~ my_scale(no.plot), data = NEON_stab_4yr.data))

all_N_4yr.result <- c()
for (i in 19:30){
  each.depend <- data_4yr[ , i]
  each.x <- data_4yr[ , 17] #no.plot
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_4yr)[i], each.result, r2)
  all_N.result <- rbind(all_N.result, d)
}
colnames(all_N_4yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_N_4yr.result, file = "all_N_effects_4yr.csv")

# 1.6 the average spatial distance of pairwise plots (spa_dist) (Supplementary Table 13)

all_spa_dis_4yr.result <- c()
for (i in 19:30){
  each.depend <- data_4yr[ , i]
  each.x <- data_4yr[ , 18] #spa_dist
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_4yr)[i], each.result, r2)
  all_spa_dis_4yr.result <- rbind(all_spa_dis_4yr.result, d)
}
colnames(all_spa_dis_4yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_spa_dis_4yr.result, file = "all_spa_dis_effects_4yr.csv")

# 1.7 the duration of data collection (duration) (Supplementary Table 14)
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div2 ~ my_scale(duration), data = NEON_stab_4yr.data))

all_duration_4yr.result <- c()
for (i in 19:30){
  each.depend <- data_4yr[ , i]
  each.x <- data_4yr[ , 16] #duration
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each_4yr.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_4yr)[i], each_4yr.result, r2)
  all_duration_4yr.result <- rbind(all_duration_4yr.result, d)
}
colnames(all_duration_4yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_duration_4yr.result, file = "all_duration_effects_4yr.csv")


#####################################################################################################
#####################################################################################################
#####################################################################################################

# Since results of the diversity-stability relationships were similar, 
# when we used the different timescales, such as >= 4-year (N=36), >= 5-year (N=24), and >= 6-year (N14). 
# Therefore, we plotted the main results in the manuscripts using the >= 4-year observation data
# in Figs. 4D-4F, Extended Data Fig. 6, Supplementary Tables 8-11. 
###############################################################################################
my_theme <- theme(legend.position = "none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(colour = "black", size = 0.5),
                  panel.border = element_blank(),
                  axis.text = element_text(size = 14, family = "Arial", color = "black"),
                  axis.title = element_text(size = 16, family = "Arial", color = "black"),
                  axis.ticks.length = unit(1, "mm"),
                  axis.line = element_line(color = "black", size = 0.5),
                  strip.text = element_text(size = 14, family = "Arial", color = "black"))
###############################################################################################

# plot 
all_effects.data <- read.csv('all_effects_4yr.csv',header=T)
variable.names(all_effects.data)

summary(all_effects.data$context)
summary(all_effects.data$factors)
summary(all_effects.data$variable)

# plot the effects on diversity (Extended Data Fig. 5)
all_div.data <- all_effects.data[grep("Biodiversity", all_effects.data$context),]
all_div.data$variable <- factor(all_div.data$variable, 
                                levels=c("tau_div","beta_div2","gamma_div","beta_div1", "alpha_div"))
all_div.data$factors <- factor(all_div.data$factors, 
                                levels=c("MAP","Precipitation seasonality","MAT","Temperature seasonality"))

effects_div.plot <- ggplot(all_div.data , aes(x = variable, y = estimate)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5, colour = "gray") +
  geom_errorbar(aes(ymin = estimate - ci, ymax = estimate + ci), size = 0.5, width= 0,alpha=1) +
  geom_point(fill="white", shape=21, size = 3.5) +
  facet_grid(.~factors) +
  coord_flip(ylim = c(-0.35,0.35)) +
  ylab(bquote(atop(paste("Effect sizes"), (NULL[paste("regression coefficients")])))) +
  xlab("Climatic factors") +
  theme_bw() +
  my_theme

# plot the effects on stability (Extended Data Fig. 5)
all_sta.data <- all_effects.data[grep("Stability", all_effects.data$context),]
all_sta.data$variable <- factor(all_sta.data$variable, 
                                levels=c("tau_sta","spa_asyn2","gamma_sta","spa_asyn1", "alpha_sta"))
all_sta.data$factors <- factor(all_sta.data$factors, 
                               levels=c("MAP","Precipitation seasonality","MAT","Temperature seasonality"))

effects_sta.plot <- ggplot(all_sta.data , aes(x = variable, y = estimate)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5, colour = "gray") +
  geom_errorbar(aes(ymin = estimate - ci, ymax = estimate + ci), size = 0.5, width= 0,alpha=1) +
  geom_point(fill="white", shape=21, size = 3.5) +
  facet_grid(.~factors) +
  coord_flip(ylim = c(-0.35,0.35)) +
  ylab(bquote(atop(paste("Effect sizes"), (NULL[paste("regression coefficients")])))) +
  xlab("Climatic factors") +
  theme_bw() +
  my_theme

# plot the effects on the diversity-stability relationship slopes at alpha, beta1, and gamma scale (Figs. 4D-4F)
all_DSRs.data <- all_effects.data[grep("DSRs", all_effects.data$context),]
all_DSRs.data$variable <- factor(all_DSRs.data$variable, 
                                levels=c("slope_gamma","slope_beta1", "slope_alpha"))
all_DSRs.data$factors <- factor(all_DSRs.data$factors, 
                               levels=c("MAP","Precipitation seasonality","MAT","Temperature seasonality"))

effects_DSRs.plot <- ggplot(all_DSRs.data , aes(x = variable, y = estimate)) +
  geom_hline(yintercept=0, linetype="dashed", size=0.5, colour = "gray") +
  geom_errorbar(aes(ymin = estimate - ci, ymax = estimate + ci), size = 0.5, width= 0,alpha=1) +
  geom_point(fill="white", shape=21, size = 3.5) +
  facet_grid(.~factors) +
  coord_flip(ylim = c(-0.35,0.35)) +
  ylab(bquote(atop(paste("Effect sizes"), (NULL[paste("regression coefficients")])))) +
  xlab("Climatic factors") +
  theme_bw() +
  my_theme


############################################
# 2 for >= 5-year observations (N=24)
NEON_stab_5yr.data <- subset(NEON_stab.data, NEON_stab.data$duration > 4)

data_5yr <- NEON_stab_5yr.data

# 2.1 mean annual precipitation (MAP) 
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(MAP_mm), data = NEON_stab_5yr.data))
summary(lm(beta_div2 ~ my_scale(MAP_mm), data = NEON_stab_5yr.data))
summary(lm(gamma_div ~ my_scale(MAP_mm), data = NEON_stab_5yr.data))
summary(lm(tau_div ~ my_scale(MAP_mm), data = NEON_stab_5yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spa_asyn1 ~ my_scale(MAP_mm), data = NEON_stab_5yr.data))
summary(lm(spa_asyn2 ~ my_scale(MAP_mm), data = NEON_stab_5yr.data))
summary(lm(tau_sta ~ my_scale(MAP_mm), data = NEON_stab_5yr.data))

all_MAP_5yr.result <- c()
for (i in 19:33){
  each.depend <- data_5yr[ , i]
  each.x <- data_5yr[ , 11] # MAP
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_5yr)[i], each.result, r2)
  all_MAP_5yr.result <- rbind(all_MAP_5yr.result, d)
}
colnames(all_MAP_5yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_MAP_5yr.result, file = "all_MAP_effects_5yr.csv")

# 2.2 mean annual temperature (MAT) 
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(MAT_C)*my_scale(MAP_mm), data = NEON_stab_5yr.data))

all_MAT_5yr.result <- c()
for (i in 19:33){
  each.depend <- data_5yr[ , i]
  each.x <- data_5yr[ , 10] # MAT
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_5yr)[i], each.result, r2)
  all_MAT_5yr.result <- rbind(all_MAT_5yr.result, d)
}
colnames(all_MAT_5yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_MAT_5yr.result, file = "all_MAT_effects_5yr.csv")

# 2.3 precipitation seasonality (Prec_s) 
# the significant (P < 0.05) effects on stability
summary(lm(spa_asyn2 ~ my_scale(Prec_s), data = NEON_stab_5yr.data))
# the significant (P < 0.05) effects on the slopes of diversity-stability relationships
summary(lm(slope_alpha_sta.alpha_div ~ my_scale(Prec_s), data = NEON_stab_5yr.data))
summary(lm(slope_gamma_sta.gamma_div ~ my_scale(Prec_s), data = NEON_stab_5yr.data))

all_Prec_s_5yr.result <- c()
for (i in 19:33){
  each.depend <- data_5yr[ , i]
  each.x <- data_5yr[ , 13] #Prec_s
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_5yr)[i], each.result, r2)
  all_Prec_s_5yr.result <- rbind(all_Prec_s_5yr.result, d)
}
colnames(all_Prec_s_5yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_Prec_s_5yr.result, file = "all_Prec_s_effects_5yr.csv")

# 2.4 temperature seasonality (Temp_s) 
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(Temp_s), data = NEON_stab_5yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spe_sta ~ my_scale(Temp_s), data = NEON_stab_5yr.data))

all_Temp_s_5yr.result <- c()
for (i in 19:33){
  each.depend <- data_5yr[ , i]
  each.x <- data_5yr[ , 12] # Temp_s
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_5yr)[i], each.result, r2)
  all_Temp_s_5yr.result <- rbind(all_Temp_s_5yr.result, d)
}
colnames(all_Temp_s_5yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_Temp_s_5yr.result, file = "all_Temp_s_effects_5yr.csv")

# 2.5 the number of plots (no.plot) 
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div2 ~ my_scale(no.plot), data = NEON_stab_5yr.data))
summary(lm(gamma_div ~ my_scale(no.plot), data = NEON_stab_5yr.data))
summary(lm(tau_div ~ my_scale(no.plot), data = NEON_stab_5yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spe_sta ~ my_scale(no.plot), data = NEON_stab_5yr.data))
# the significant (P < 0.05) effects on the slopes of diversity-stability relationships
summary(lm(slope_alpha_sta.alpha_div ~ my_scale(no.plot), data = NEON_stab_5yr.data))

all_N_5yr.result <- c()
for (i in 19:33){
  each.depend <- data_5yr[ , i]
  each.x <- data_5yr[ , 17] #no.plot
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_5yr)[i], each.result, r2)
  all_N_5yr.result <- rbind(all_N_5yr.result, d)
}
colnames(all_N_5yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_N_5yr.result, file = "all_N_effects_5yr.csv")

# 2.6 the average spatial distance of pairwise plots (spa_dist) 

all_spa_dis_5yr.result <- c()
for (i in 19:33){
  each.depend <- data_5yr[ , i]
  each.x <- data_5yr[ , 18] #spa_dist
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_5yr)[i], each.result, r2)
  all_spa_dis_5yr.result <- rbind(all_spa_dis_5yr.result, d)
}
colnames(all_spa_dis_5yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_spa_dis_5yr.result, file = "all_spa_dis_effects_5yr.csv")

# 2.7 the duration of data collection (duration) 

all_duration_5yr.result <- c()
for (i in 19:33){
  each.depend <- data_5yr[ , i]
  each.x <- data_5yr[ , 16] #duration
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each_5yr.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_5yr)[i], each_5yr.result, r2)
  all_duration_5yr.result <- rbind(all_duration_5yr.result, d)
}
colnames(all_duration_5yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_duration_5yr.result, file = "all_duration_effects_5yr.csv")



############################################
# 3 for >= 6-year observations (N=14)
NEON_stab_6yr.data <- subset(NEON_stab.data, NEON_stab.data$duration > 5)

data_6yr <- NEON_stab_6yr.data

# 3.1 mean annual precipitation (MAP) 
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(MAP_mm), data = NEON_stab_6yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spa_asyn1 ~ my_scale(MAP_mm), data = NEON_stab_6yr.data))

all_MAP_6yr.result <- c()
for (i in 19:33){
  each.depend <- data_6yr[ , i]
  each.x <- data_6yr[ , 11] # MAP
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_6yr)[i], each.result, r2)
  all_MAP_6yr.result <- rbind(all_MAP_6yr.result, d)
}
colnames(all_MAP_6yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_MAP_6yr.result, file = "all_MAP_effects_6yr.csv")

# 3.2 mean annual temperature (MAT) 

all_MAT_6yr.result <- c()
for (i in 19:33){
  each.depend <- data_6yr[ , i]
  each.x <- data_6yr[ , 10] # MAT
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_6yr)[i], each.result, r2)
  all_MAT_6yr.result <- rbind(all_MAT_6yr.result, d)
}
colnames(all_MAT_6yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_MAT_6yr.result, file = "all_MAT_effects_6yr.csv")

# 3.3 precipitation seasonality (Prec_s) 
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div1 ~ my_scale(Prec_s), data = NEON_stab_6yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spa_asyn1 ~ my_scale(Prec_s), data = NEON_stab_6yr.data))

all_Prec_s_6yr.result <- c()
for (i in 19:33){
  each.depend <- data_6yr[ , i]
  each.x <- data_6yr[ , 13] #Prec_s
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_6yr)[i], each.result, r2)
  all_Prec_s_6yr.result <- rbind(all_Prec_s_6yr.result, d)
}
colnames(all_Prec_s_6yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_Prec_s_6yr.result, file = "all_Prec_s_effects_6yr.csv")

# 3.4 temperature seasonality (Temp_s) 

all_Temp_s_6yr.result <- c()
for (i in 19:33){
  each.depend <- data_6yr[ , i]
  each.x <- data_6yr[ , 12] # Temp_s
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_6yr)[i], each.result, r2)
  all_Temp_s_6yr.result <- rbind(all_Temp_s_6yr.result, d)
}
colnames(all_Temp_s_6yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_Temp_s_6yr.result, file = "all_Temp_s_effects_6yr.csv")

# 3.5 the number of plots (no.plot) 
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div2 ~ my_scale(no.plot), data = NEON_stab_6yr.data))
# the significant (P < 0.05) effects on stability
summary(lm(spe_sta ~ my_scale(no.plot), data = NEON_stab_6yr.data))

all_N_6yr.result <- c()
for (i in 19:33){
  each.depend <- data_6yr[ , i]
  each.x <- data_6yr[ , 17] #no.plot
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_6yr)[i], each.result, r2)
  all_N_6yr.result <- rbind(all_N_6yr.result, d)
}
colnames(all_N_6yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_N_6yr.result, file = "all_N_effects_6yr.csv")

# 3.6 the average spatial distance of pairwise plots (spa_dist) 

all_spa_dis_6yr.result <- c()
for (i in 19:33){
  each.depend <- data_6yr[ , i]
  each.x <- data_6yr[ , 18] #spa_dist
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_6yr)[i], each.result, r2)
  all_spa_dis_6yr.result <- rbind(all_spa_dis_6yr.result, d)
}
colnames(all_spa_dis_6yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_spa_dis_6yr.result, file = "all_spa_dis_effects_6yr.csv")

# 3.7 the duration of data collection (duration) 
# the significant (P < 0.05) effects on diversity
summary(lm(beta_div2 ~ my_scale(duration), data = NEON_stab_6yr.data))

all_duration_6yr.result <- c()
for (i in 19:33){
  each.depend <- data_6yr[ , i]
  each.x <- data_6yr[ , 16] #duration
  each.fit <- lm(each.depend ~ my_scale(each.x))
  each_6yr.result <- summary(each.fit)$coefficients[2, 1:4]; r2 <- r.squaredGLMM(each.fit)
  d <- c(colnames(data_6yr)[i], each_6yr.result, r2)
  all_duration_6yr.result <- rbind(all_duration_6yr.result, d)
}
colnames(all_duration_6yr.result) <- c("variable", "estimate", "se", "t value", "p", "r2m", "r2c")
#write.csv(all_duration_6yr.result, file = "all_duration_effects_6yr.csv")
