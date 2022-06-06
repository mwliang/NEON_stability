# Load packages
library(nlme)
library(lme4)
library(MuMIn)
library(ggplot2)

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
##############################################################################################################

############################################
#transformation function
scaleFUN <- function(x) sprintf("%.1f", x)
############################################

###############################################################################################
my_theme <- theme(legend.position = "none",
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(colour = "black", size = 0.5),
                  panel.border = element_rect(colour = "black", size = 0.5),
                  axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
                  axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
                  axis.title.y = element_text(size = 18, family = "Arial", color = "black"),
                  axis.title.x = element_text(size = 18, family = "Arial", color = "black"),
                  axis.ticks.length = unit(1, "mm"),
                  axis.line = element_line(color = "black", size = 0.5),
                  plot.tag = element_text(size = 20, family = "Arial", color = "black")) 
###############################################################################################


############################################
# 1 for >= 4-year observations (N=945)
NEON_stab_plots_4yr.data <- NEON_stab_plots.data

# 1.1 the diversity-stability relationships at alpha scale 
# fit the model examining the effects of alpha diversity on alpha stability (Fig. 4A and Supplementary Table 3)
alpha_alpha_sta_4yr_lmm.fit <- lme(alpha_sta ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                   correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_4yr.data)
summary(alpha_alpha_sta_4yr_lmm.fit)
anova(alpha_alpha_sta_4yr_lmm.fit)
r.squaredGLMM(alpha_alpha_sta_4yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_alpha_sta_4yr_lmm.fit), NEON_stab_plots_4yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_alpha_sta_4yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_alpha_sta_4yr_lmm.fit, which=1:4)

############################################## generating new data ###############################################
alpha_alpha_sta_4yr_newdat.lme = data.frame(siteID = NEON_stab_plots_4yr.data$siteID,
                                            alpha_div = NEON_stab_plots_4yr.data$alpha_div)
head(alpha_alpha_sta_4yr_newdat.lme)
alpha_alpha_sta_4yr_newdat.lme$predlme = predict(alpha_alpha_sta_4yr_lmm.fit, newdata = alpha_alpha_sta_4yr_newdat.lme, level = 0)
alpha_alpha_sta_4yr_des = model.matrix(formula(alpha_alpha_sta_4yr_lmm.fit)[-2], alpha_alpha_sta_4yr_newdat.lme)
alpha_alpha_sta_4yr_predvar = diag( alpha_alpha_sta_4yr_des %*% vcov(alpha_alpha_sta_4yr_lmm.fit) %*% t(alpha_alpha_sta_4yr_des))
alpha_alpha_sta_4yr_newdat.lme$lower = with(alpha_alpha_sta_4yr_newdat.lme, predlme - 2*sqrt(alpha_alpha_sta_4yr_predvar))
alpha_alpha_sta_4yr_newdat.lme$upper = with(alpha_alpha_sta_4yr_newdat.lme, predlme + 2*sqrt(alpha_alpha_sta_4yr_predvar))
################################################# plotting ######################################################
alpha_alpha_sta_4yr.plot <- ggplot(NEON_stab_plots_4yr.data, aes(y = alpha_sta, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = alpha_alpha_sta_4yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_alpha_sta_4yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  annotate("text", x = 2, y = -0.1, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.10 ***\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 2, y = -0.45, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.36 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 1.2 the diversity-stability relationships at beta1 scale 
# fit the model examining the effects of beta diversity 1 on spatial asynchrony 1 (beta stability 1) (Fig. 4B and Supplementary Table 3)
beta1_spa_asy1_4yr_lmm.fit <- lme(spa_asy1 ~ beta_div1, random = ~1|siteID/nlcdClass, 
                                  correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_4yr.data)
summary(beta1_spa_asy1_4yr_lmm.fit)
anova(beta1_spa_asy1_4yr_lmm.fit)
r.squaredGLMM(beta1_spa_asy1_4yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(beta1_spa_asy1_4yr_lmm.fit), NEON_stab_plots_4yr.data$beta_div1)
# Assess normality of residuals
hist(beta1_spa_asy1_4yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(beta1_spa_asy1_4yr_lmm.fit)

############################################## generating new data ###############################################
beta1_spa_asy1_4yr_newdat.lme = data.frame(siteID = NEON_stab_plots_4yr.data$siteID,
                                           beta_div1 = NEON_stab_plots_4yr.data$beta_div1)
head(beta1_spa_asy1_4yr_newdat.lme)
beta1_spa_asy1_4yr_newdat.lme$predlme = predict(beta1_spa_asy1_4yr_lmm.fit, newdata = beta1_spa_asy1_4yr_newdat.lme, level = 0)
beta1_spa_asy1_4yr_des = model.matrix(formula(beta1_spa_asy1_4yr_lmm.fit)[-2], beta1_spa_asy1_4yr_newdat.lme)
beta1_spa_asy1_4yr_predvar = diag( beta1_spa_asy1_4yr_des %*% vcov(beta1_spa_asy1_4yr_lmm.fit) %*% t(beta1_spa_asy1_4yr_des))
beta1_spa_asy1_4yr_newdat.lme$lower = with(beta1_spa_asy1_4yr_newdat.lme, predlme - 2*sqrt(beta1_spa_asy1_4yr_predvar))
beta1_spa_asy1_4yr_newdat.lme$upper = with(beta1_spa_asy1_4yr_newdat.lme, predlme + 2*sqrt(beta1_spa_asy1_4yr_predvar))
################################################# plotting ######################################################
beta1_spa_asy1_4yr.plot <- ggplot(NEON_stab_plots_4yr.data, aes(y = spa_asy1, x = beta_div1, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = beta1_spa_asy1_4yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), 
              alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = beta1_spa_asy1_4yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  #xlim(0, 2) +
  annotate("text", x = 0.2, y = 3.0, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.02 **\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 0.2, y = 2.5, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.25 **\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 1.3 the diversity-stability relationships at gamma scale 
# fit the model examining the effects of gamma diversity on gamma stability (Fig. 4C and Supplementary Table 3)
gamma_gamma_sta_4yr_lmm.fit <- lme(gamma_sta ~ gamma_div, random = ~1|siteID/nlcdClass, 
                                   correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_4yr.data)
summary(gamma_gamma_sta_4yr_lmm.fit)
anova(gamma_gamma_sta_4yr_lmm.fit)
r.squaredGLMM(gamma_gamma_sta_4yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(gamma_gamma_sta_4yr_lmm.fit), NEON_stab_plots_4yr.data$gamma_div)
# Assess normality of residuals
hist(gamma_gamma_sta_4yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(gamma_gamma_sta_4yr_lmm.fit)

############################################## generating new data ###############################################
gamma_gamma_sta_4yr_newdat.lme = data.frame(siteID = NEON_stab_plots_4yr.data$siteID,
                                            gamma_div = NEON_stab_plots_4yr.data$gamma_div)
head(gamma_gamma_sta_4yr_newdat.lme)
gamma_gamma_sta_4yr_newdat.lme$predlme = predict(gamma_gamma_sta_4yr_lmm.fit, newdata = gamma_gamma_sta_4yr_newdat.lme, level = 0)
gamma_gamma_sta_4yr_des = model.matrix(formula(gamma_gamma_sta_4yr_lmm.fit)[-2], gamma_gamma_sta_4yr_newdat.lme)
gamma_gamma_sta_4yr_predvar = diag( gamma_gamma_sta_4yr_des %*% vcov(gamma_gamma_sta_4yr_lmm.fit) %*% t(gamma_gamma_sta_4yr_des))
gamma_gamma_sta_4yr_newdat.lme$lower = with(gamma_gamma_sta_4yr_newdat.lme, predlme - 2*sqrt(gamma_gamma_sta_4yr_predvar))
gamma_gamma_sta_4yr_newdat.lme$upper = with(gamma_gamma_sta_4yr_newdat.lme, predlme + 2*sqrt(gamma_gamma_sta_4yr_predvar))
################################################# plotting ######################################################
gamma_gamma_sta_4yr.plot <- ggplot(NEON_stab_plots_4yr.data, aes(y = gamma_sta, x = gamma_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = gamma_gamma_sta_4yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), 
              alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = gamma_gamma_sta_4yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  annotate("text", x = 0.9, y = 4.3, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.03 ***\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 0.9, y = 3.5, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.31 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 1.4 the diversity-stability relationships at population-level 
# fit the model examining the effects of alpha diversity on species stability (Supplementary Fig. 1A and Table 4)
alpha_spe_sta_4yr_lmm.fit <- lme(spe_sta ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                 correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_4yr.data)
summary(alpha_spe_sta_4yr_lmm.fit)
anova(alpha_spe_sta_4yr_lmm.fit)
r.squaredGLMM(alpha_spe_sta_4yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_spe_sta_4yr_lmm.fit), NEON_stab_plots_4yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_spe_sta_4yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_spe_sta_4yr_lmm.fit)

############################################## generating new data ###############################################
alpha_spe_sta_4yr_newdat.lme = data.frame(siteID = NEON_stab_plots_4yr.data$siteID,
                                          alpha_div = NEON_stab_plots_4yr.data$alpha_div)
head(alpha_spe_sta_4yr_newdat.lme)
alpha_spe_sta_4yr_newdat.lme$predlme = predict(alpha_spe_sta_4yr_lmm.fit, newdata = alpha_spe_sta_4yr_newdat.lme, level = 0)
alpha_spe_sta_4yr_des = model.matrix(formula(alpha_spe_sta_4yr_lmm.fit)[-2], alpha_spe_sta_4yr_newdat.lme)
alpha_spe_sta_4yr_predvar = diag( alpha_spe_sta_4yr_des %*% vcov(alpha_spe_sta_4yr_lmm.fit) %*% t(alpha_spe_sta_4yr_des))
alpha_spe_sta_4yr_newdat.lme$lower = with(alpha_spe_sta_4yr_newdat.lme, predlme - 2*sqrt(alpha_spe_sta_4yr_predvar))
alpha_spe_sta_4yr_newdat.lme$upper = with(alpha_spe_sta_4yr_newdat.lme, predlme + 2*sqrt(alpha_spe_sta_4yr_predvar))
################################################# plotting ######################################################
alpha_spe_sta_4yr.plot <- ggplot(NEON_stab_plots_4yr.data, aes(y = spe_sta, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  #geom_ribbon(data = alpha_spa_sta_4yr_newdat.lme, 
  #aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  #scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_spe_sta_4yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black", linetype ="dashed") +
  ylab(bquote(paste("Species stability ("*NULL[paste(log[10])]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  annotate("text", x = 2.1, y = 1.0, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.001\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 2.1, y = 0.75, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.53\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 1.5 the diversity-stability relationships at population-level 
# fit the model examining the effects of alpha diversity on species asynchrony (Supplementary Fig. 1B and Table 4)
alpha_spe_asy_4yr_lmm.fit <- lme(spe_asy ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                 correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_4yr.data)
summary(alpha_spe_asy_4yr_lmm.fit)
anova(alpha_spe_asy_4yr_lmm.fit)
r.squaredGLMM(alpha_spe_asy_4yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_spe_asy_4yr_lmm.fit), NEON_stab_plots_4yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_spe_asy_4yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_spe_asy_4yr_lmm.fit)

############################################## generating new data ###############################################
alpha_spe_asy_4yr_newdat.lme = data.frame(siteID = NEON_stab_plots_4yr.data$siteID,
                                          alpha_div = NEON_stab_plots_4yr.data$alpha_div)
head(alpha_spe_asy_4yr_newdat.lme)
alpha_spe_asy_4yr_newdat.lme$predlme = predict(alpha_spe_asy_4yr_lmm.fit, newdata = alpha_spe_asy_4yr_newdat.lme, level = 0)
alpha_spe_asy_4yr_des = model.matrix(formula(alpha_spe_asy_4yr_lmm.fit)[-2], alpha_spe_asy_4yr_newdat.lme)
alpha_spe_asy_4yr_predvar = diag( alpha_spe_asy_4yr_des %*% vcov(alpha_spe_asy_4yr_lmm.fit) %*% t(alpha_spe_asy_4yr_des))
alpha_spe_asy_4yr_newdat.lme$lower = with(alpha_spe_asy_4yr_newdat.lme, predlme - 2*sqrt(alpha_spe_asy_4yr_predvar))
alpha_spe_asy_4yr_newdat.lme$upper = with(alpha_spe_asy_4yr_newdat.lme, predlme + 2*sqrt(alpha_spe_asy_4yr_predvar))
################################################# plotting ######################################################
alpha_spe_asy_4yr.plot <- ggplot(NEON_stab_plots_4yr.data, aes(y = spe_asy, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = alpha_spe_asy_4yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_spe_asy_4yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste("Species asynchrony ("*NULL[paste(log[10])]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  annotate("text", x = 1.5, y = 2.1, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.10 ***\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 1.5, y = 1.8, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.59 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme


############################################
# 2 for >= 5-year observations (N=570)
NEON_stab_plots_5yr.data <- subset(NEON_stab_plots.data, NEON_stab_plots.data$duration > 4)

# 2.1 the diversity-stability relationships at alpha scale 
# fit the model examining the effects of alpha diversity on alpha stability  (Supplementary Table 3)
alpha_alpha_sta_5yr_lmm.fit <- lme(alpha_sta ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                   correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_5yr.data)
summary(alpha_alpha_sta_5yr_lmm.fit)
anova(alpha_alpha_sta_5yr_lmm.fit)
r.squaredGLMM(alpha_alpha_sta_5yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_alpha_sta_5yr_lmm.fit), NEON_stab_plots_5yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_alpha_sta_5yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_alpha_sta_5yr_lmm.fit, which=1:4)

############################################## generating new data ###############################################
alpha_alpha_sta_5yr_newdat.lme = data.frame(siteID = NEON_stab_plots_5yr.data$siteID,
                                            alpha_div = NEON_stab_plots_5yr.data$alpha_div)
head(alpha_alpha_sta_5yr_newdat.lme)
alpha_alpha_sta_5yr_newdat.lme$predlme = predict(alpha_alpha_sta_5yr_lmm.fit, newdata = alpha_alpha_sta_5yr_newdat.lme, level = 0)
alpha_alpha_sta_5yr_des = model.matrix(formula(alpha_alpha_sta_5yr_lmm.fit)[-2], alpha_alpha_sta_5yr_newdat.lme)
alpha_alpha_sta_5yr_predvar = diag( alpha_alpha_sta_5yr_des %*% vcov(alpha_alpha_sta_5yr_lmm.fit) %*% t(alpha_alpha_sta_5yr_des))
alpha_alpha_sta_5yr_newdat.lme$lower = with(alpha_alpha_sta_5yr_newdat.lme, predlme - 2*sqrt(alpha_alpha_sta_5yr_predvar))
alpha_alpha_sta_5yr_newdat.lme$upper = with(alpha_alpha_sta_5yr_newdat.lme, predlme + 2*sqrt(alpha_alpha_sta_5yr_predvar))
################################################# plotting ######################################################
alpha_alpha_sta_5yr.plot <- ggplot(NEON_stab_plots_5yr.data, aes(y = alpha_sta, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = alpha_alpha_sta_5yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_alpha_sta_5yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  annotate("text", x = 2, y = -0.1, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.13 ***\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 2, y = -0.45, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.37 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 2.2 the diversity-stability relationships at beta1 scale 
# fit the model examining the effects of beta diversity 1 on spatial asynchrony 1  (Supplementary Table 3)
beta1_spa_asy1_5yr_lmm.fit <- lme(spa_asy1 ~ beta_div1, random = ~1|siteID/nlcdClass, 
                                  correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_5yr.data)
summary(beta1_spa_asy1_5yr_lmm.fit)
anova(beta1_spa_asy1_5yr_lmm.fit)
r.squaredGLMM(beta1_spa_asy1_5yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(beta1_spa_asy1_5yr_lmm.fit), NEON_stab_plots_5yr.data$beta_div1)
# Assess normality of residuals
hist(beta1_spa_asy1_5yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(beta1_spa_asy1_5yr_lmm.fit)

############################################## generating new data ###############################################
beta1_spa_asy1_5yr_newdat.lme = data.frame(siteID = NEON_stab_plots_5yr.data$siteID,
                                           beta_div1 = NEON_stab_plots_5yr.data$beta_div1)
head(beta1_spa_asy1_5yr_newdat.lme)
beta1_spa_asy1_5yr_newdat.lme$predlme = predict(beta1_spa_asy1_5yr_lmm.fit, newdata = beta1_spa_asy1_5yr_newdat.lme, level = 0)
beta1_spa_asy1_5yr_des = model.matrix(formula(beta1_spa_asy1_5yr_lmm.fit)[-2], beta1_spa_asy1_5yr_newdat.lme)
beta1_spa_asy1_5yr_predvar = diag( beta1_spa_asy1_5yr_des %*% vcov(beta1_spa_asy1_5yr_lmm.fit) %*% t(beta1_spa_asy1_5yr_des))
beta1_spa_asy1_5yr_newdat.lme$lower = with(beta1_spa_asy1_5yr_newdat.lme, predlme - 2*sqrt(beta1_spa_asy1_5yr_predvar))
beta1_spa_asy1_5yr_newdat.lme$upper = with(beta1_spa_asy1_5yr_newdat.lme, predlme + 2*sqrt(beta1_spa_asy1_5yr_predvar))
################################################# plotting ######################################################
beta1_spa_asy1_5yr.plot <- ggplot(NEON_stab_plots_5yr.data, aes(y = spa_asy1, x = beta_div1, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = beta1_spa_asy1_5yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), 
              alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = beta1_spa_asy1_5yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  #xlim(0, 2) +
  annotate("text", x = 0.2, y = 3.0, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.01 *\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 0.2, y = 2.5, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.27 *\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 2.3 the diversity-stability relationships at gamma scale 
# fit the model examining the effects of gamma diversity on gamma stability (Supplementary Table 3)
gamma_gamma_sta_5yr_lmm.fit <- lme(gamma_sta ~ gamma_div, random = ~1|siteID/nlcdClass, 
                                   correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_5yr.data)
summary(gamma_gamma_sta_5yr_lmm.fit)
anova(gamma_gamma_sta_5yr_lmm.fit)
r.squaredGLMM(gamma_gamma_sta_5yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(gamma_gamma_sta_5yr_lmm.fit), NEON_stab_plots_5yr.data$gamma_div)
# Assess normality of residuals
hist(gamma_gamma_sta_5yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(gamma_gamma_sta_5yr_lmm.fit)

############################################## generating new data ###############################################
gamma_gamma_sta_5yr_newdat.lme = data.frame(siteID = NEON_stab_plots_5yr.data$siteID,
                                            gamma_div = NEON_stab_plots_5yr.data$gamma_div)
head(gamma_gamma_sta_5yr_newdat.lme)
gamma_gamma_sta_5yr_newdat.lme$predlme = predict(gamma_gamma_sta_5yr_lmm.fit, newdata = gamma_gamma_sta_5yr_newdat.lme, level = 0)
gamma_gamma_sta_5yr_des = model.matrix(formula(gamma_gamma_sta_5yr_lmm.fit)[-2], gamma_gamma_sta_5yr_newdat.lme)
gamma_gamma_sta_5yr_predvar = diag( gamma_gamma_sta_5yr_des %*% vcov(gamma_gamma_sta_5yr_lmm.fit) %*% t(gamma_gamma_sta_5yr_des))
gamma_gamma_sta_5yr_newdat.lme$lower = with(gamma_gamma_sta_5yr_newdat.lme, predlme - 2*sqrt(gamma_gamma_sta_5yr_predvar))
gamma_gamma_sta_5yr_newdat.lme$upper = with(gamma_gamma_sta_5yr_newdat.lme, predlme + 2*sqrt(gamma_gamma_sta_5yr_predvar))
################################################# plotting ######################################################
gamma_gamma_sta_5yr.plot <- ggplot(NEON_stab_plots_5yr.data, aes(y = gamma_sta, x = gamma_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = gamma_gamma_sta_5yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), 
              alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = gamma_gamma_sta_5yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  annotate("text", x = 0.9, y = 4.3, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.05 ***\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 0.9, y = 3.5, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.34 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 2.4 the diversity-stability relationships at population-level 
# fit the model examining the effects of alpha diversity on species stability 
alpha_spe_sta_5yr_lmm.fit <- lme(spe_sta ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                 correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_5yr.data)
summary(alpha_spe_sta_5yr_lmm.fit)
anova(alpha_spe_sta_5yr_lmm.fit)
r.squaredGLMM(alpha_spe_sta_5yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_spe_sta_5yr_lmm.fit), NEON_stab_plots_5yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_spe_sta_5yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_spe_sta_5yr_lmm.fit)

############################################## generating new data ###############################################
alpha_spe_sta_5yr_newdat.lme = data.frame(siteID = NEON_stab_plots_5yr.data$siteID,
                                          alpha_div = NEON_stab_plots_5yr.data$alpha_div)
head(alpha_spe_sta_5yr_newdat.lme)
alpha_spe_sta_5yr_newdat.lme$predlme = predict(alpha_spe_sta_5yr_lmm.fit, newdata = alpha_spe_sta_5yr_newdat.lme, level = 0)
alpha_spe_sta_5yr_des = model.matrix(formula(alpha_spe_sta_5yr_lmm.fit)[-2], alpha_spe_sta_5yr_newdat.lme)
alpha_spe_sta_5yr_predvar = diag( alpha_spe_sta_5yr_des %*% vcov(alpha_spe_sta_5yr_lmm.fit) %*% t(alpha_spe_sta_5yr_des))
alpha_spe_sta_5yr_newdat.lme$lower = with(alpha_spe_sta_5yr_newdat.lme, predlme - 2*sqrt(alpha_spe_sta_5yr_predvar))
alpha_spe_sta_5yr_newdat.lme$upper = with(alpha_spe_sta_5yr_newdat.lme, predlme + 2*sqrt(alpha_spe_sta_5yr_predvar))
################################################# plotting ######################################################
alpha_spe_sta_5yr.plot <- ggplot(NEON_stab_plots_5yr.data, aes(y = spe_sta, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  #geom_ribbon(data = alpha_spa_sta_5yr_newdat.lme, 
  #aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  #scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_spe_sta_5yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black", linetype ="dashed") +
  ylab(bquote(paste("Species stability ("*NULL[paste(log[10])]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "D") +
  annotate("text", x = 2.1, y = 1.0, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.004\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 2.1, y = 0.75, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.43\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 2.5 the diversity-stability relationships at population-level 
# fit the model examining the effects of alpha diversity on species asynchrony 
alpha_spe_asy_5yr_lmm.fit <- lme(spe_asy ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                 correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_5yr.data)
summary(alpha_spe_asy_5yr_lmm.fit)
anova(alpha_spe_asy_5yr_lmm.fit)
r.squaredGLMM(alpha_spe_asy_5yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_spe_asy_5yr_lmm.fit), NEON_stab_plots_5yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_spe_asy_5yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_spe_asy_5yr_lmm.fit)

############################################## generating new data ###############################################
alpha_spe_asy_5yr_newdat.lme = data.frame(siteID = NEON_stab_plots_5yr.data$siteID,
                                          alpha_div = NEON_stab_plots_5yr.data$alpha_div)
head(alpha_spe_asy_5yr_newdat.lme)
alpha_spe_asy_5yr_newdat.lme$predlme = predict(alpha_spe_asy_5yr_lmm.fit, newdata = alpha_spe_asy_5yr_newdat.lme, level = 0)
alpha_spe_asy_5yr_des = model.matrix(formula(alpha_spe_asy_5yr_lmm.fit)[-2], alpha_spe_asy_5yr_newdat.lme)
alpha_spe_asy_5yr_predvar = diag( alpha_spe_asy_5yr_des %*% vcov(alpha_spe_asy_5yr_lmm.fit) %*% t(alpha_spe_asy_5yr_des))
alpha_spe_asy_5yr_newdat.lme$lower = with(alpha_spe_asy_5yr_newdat.lme, predlme - 2*sqrt(alpha_spe_asy_5yr_predvar))
alpha_spe_asy_5yr_newdat.lme$upper = with(alpha_spe_asy_5yr_newdat.lme, predlme + 2*sqrt(alpha_spe_asy_5yr_predvar))
################################################# plotting ######################################################
alpha_spe_asy_5yr.plot <- ggplot(NEON_stab_plots_5yr.data, aes(y = spe_asy, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = alpha_spe_asy_5yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_spe_asy_5yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste("Species asynchrony ("*NULL[paste(log[10])]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "E") +
  annotate("text", x = 1.5, y = 2.1, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.14 ***\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 1.5, y = 1.8, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.60 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

############################################
# 3 for >= 6-year observations (N=309)
NEON_stab_plots_6yr.data <- subset(NEON_stab_plots.data, NEON_stab_plots.data$duration > 5)

# 3.1 the diversity-stability relationships at alpha scale 
# fit the model examining the effects of alpha diversity on alpha stability (Supplementary Table 3)
alpha_alpha_sta_6yr_lmm.fit <- lme(alpha_sta ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                   correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_6yr.data)
summary(alpha_alpha_sta_6yr_lmm.fit)
anova(alpha_alpha_sta_6yr_lmm.fit)
r.squaredGLMM(alpha_alpha_sta_6yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_alpha_sta_6yr_lmm.fit), NEON_stab_plots_6yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_alpha_sta_6yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_alpha_sta_6yr_lmm.fit, which=1:4)

############################################## generating new data ###############################################
alpha_alpha_sta_6yr_newdat.lme = data.frame(siteID = NEON_stab_plots_6yr.data$siteID,
                                            alpha_div = NEON_stab_plots_6yr.data$alpha_div)
head(alpha_alpha_sta_6yr_newdat.lme)
alpha_alpha_sta_6yr_newdat.lme$predlme = predict(alpha_alpha_sta_6yr_lmm.fit, newdata = alpha_alpha_sta_6yr_newdat.lme, level = 0)
alpha_alpha_sta_6yr_des = model.matrix(formula(alpha_alpha_sta_6yr_lmm.fit)[-2], alpha_alpha_sta_6yr_newdat.lme)
alpha_alpha_sta_6yr_predvar = diag( alpha_alpha_sta_6yr_des %*% vcov(alpha_alpha_sta_6yr_lmm.fit) %*% t(alpha_alpha_sta_6yr_des))
alpha_alpha_sta_6yr_newdat.lme$lower = with(alpha_alpha_sta_6yr_newdat.lme, predlme - 2*sqrt(alpha_alpha_sta_6yr_predvar))
alpha_alpha_sta_6yr_newdat.lme$upper = with(alpha_alpha_sta_6yr_newdat.lme, predlme + 2*sqrt(alpha_alpha_sta_6yr_predvar))
################################################# plotting ######################################################
alpha_alpha_sta_6yr.plot <- ggplot(NEON_stab_plots_6yr.data, aes(y = alpha_sta, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = alpha_alpha_sta_6yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_alpha_sta_6yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  annotate("text", x = 2, y = -0.1, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.15 ***\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 2, y = -0.45, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.41 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 3.2 the diversity-stability relationships at beta1 scale 
# fit the model examining the effects of beta diversity 1 on spatial asynchrony 1 (Supplementary Table 3)
beta1_spa_asy1_6yr_lmm.fit <- lme(spa_asy1 ~ beta_div1, random = ~1|siteID/nlcdClass, 
                                  correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_6yr.data)
summary(beta1_spa_asy1_6yr_lmm.fit)
anova(beta1_spa_asy1_6yr_lmm.fit)
r.squaredGLMM(beta1_spa_asy1_6yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(beta1_spa_asy1_6yr_lmm.fit), NEON_stab_plots_6yr.data$beta_div1)
# Assess normality of residuals
hist(beta1_spa_asy1_6yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(beta1_spa_asy1_6yr_lmm.fit)

############################################## generating new data ###############################################
beta1_spa_asy1_6yr_newdat.lme = data.frame(siteID = NEON_stab_plots_6yr.data$siteID,
                                           beta_div1 = NEON_stab_plots_6yr.data$beta_div1)
head(beta1_spa_asy1_6yr_newdat.lme)
beta1_spa_asy1_6yr_newdat.lme$predlme = predict(beta1_spa_asy1_6yr_lmm.fit, newdata = beta1_spa_asy1_6yr_newdat.lme, level = 0)
beta1_spa_asy1_6yr_des = model.matrix(formula(beta1_spa_asy1_6yr_lmm.fit)[-2], beta1_spa_asy1_6yr_newdat.lme)
beta1_spa_asy1_6yr_predvar = diag( beta1_spa_asy1_6yr_des %*% vcov(beta1_spa_asy1_6yr_lmm.fit) %*% t(beta1_spa_asy1_6yr_des))
beta1_spa_asy1_6yr_newdat.lme$lower = with(beta1_spa_asy1_6yr_newdat.lme, predlme - 2*sqrt(beta1_spa_asy1_6yr_predvar))
beta1_spa_asy1_6yr_newdat.lme$upper = with(beta1_spa_asy1_6yr_newdat.lme, predlme + 2*sqrt(beta1_spa_asy1_6yr_predvar))
################################################# plotting ######################################################
beta1_spa_asy1_6yr.plot <- ggplot(NEON_stab_plots_6yr.data, aes(y = spa_asy1, x = beta_div1, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  #geom_ribbon(data = beta1_spa_asy1_6yr_newdat.lme, 
  #aes(y = NULL, ymin = lower, ymax = upper, color = NULL), 
  #alpha = .15) +  # add Confidence interval
  #scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = beta1_spa_asy1_6yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black", linetype ="dashed") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  #xlim(0, 2) +
  annotate("text", x = 0.2, y = 3.0, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.009\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 0.2, y = 2.5, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.33\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 3.3 the diversity-stability relationships at gamma scale 
# fit the model examining the effects of gamma diversity on gamma stability (Supplementary Table 3)
gamma_gamma_sta_6yr_lmm.fit <- lme(gamma_sta ~ gamma_div, random = ~1|siteID/nlcdClass, 
                                   correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_6yr.data)
summary(gamma_gamma_sta_6yr_lmm.fit)
anova(gamma_gamma_sta_6yr_lmm.fit)
r.squaredGLMM(gamma_gamma_sta_6yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(gamma_gamma_sta_6yr_lmm.fit), NEON_stab_plots_6yr.data$gamma_div)
# Assess normality of residuals
hist(gamma_gamma_sta_6yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(gamma_gamma_sta_6yr_lmm.fit)

############################################## generating new data ###############################################
gamma_gamma_sta_6yr_newdat.lme = data.frame(siteID = NEON_stab_plots_6yr.data$siteID,
                                            gamma_div = NEON_stab_plots_6yr.data$gamma_div)
head(gamma_gamma_sta_6yr_newdat.lme)
gamma_gamma_sta_6yr_newdat.lme$predlme = predict(gamma_gamma_sta_6yr_lmm.fit, newdata = gamma_gamma_sta_6yr_newdat.lme, level = 0)
gamma_gamma_sta_6yr_des = model.matrix(formula(gamma_gamma_sta_6yr_lmm.fit)[-2], gamma_gamma_sta_6yr_newdat.lme)
gamma_gamma_sta_6yr_predvar = diag( gamma_gamma_sta_6yr_des %*% vcov(gamma_gamma_sta_6yr_lmm.fit) %*% t(gamma_gamma_sta_6yr_des))
gamma_gamma_sta_6yr_newdat.lme$lower = with(gamma_gamma_sta_6yr_newdat.lme, predlme - 2*sqrt(gamma_gamma_sta_6yr_predvar))
gamma_gamma_sta_6yr_newdat.lme$upper = with(gamma_gamma_sta_6yr_newdat.lme, predlme + 2*sqrt(gamma_gamma_sta_6yr_predvar))
################################################# plotting ######################################################
gamma_gamma_sta_6yr.plot <- ggplot(NEON_stab_plots_6yr.data, aes(y = gamma_sta, x = gamma_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = gamma_gamma_sta_6yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), 
              alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = gamma_gamma_sta_6yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  annotate("text", x = 0.9, y = 4.3, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.04 **\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 0.9, y = 3.5, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.36 **\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 3.4 the diversity-stability relationships at population-level 
# fit the model examining the effects of alpha diversity on species stability 
alpha_spe_sta_6yr_lmm.fit <- lme(spe_sta ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                 correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_6yr.data)
summary(alpha_spe_sta_6yr_lmm.fit)
anova(alpha_spe_sta_6yr_lmm.fit)
r.squaredGLMM(alpha_spe_sta_6yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_spe_sta_6yr_lmm.fit), NEON_stab_plots_6yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_spe_sta_6yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_spe_sta_6yr_lmm.fit)

############################################## generating new data ###############################################
alpha_spe_sta_6yr_newdat.lme = data.frame(siteID = NEON_stab_plots_6yr.data$siteID,
                                          alpha_div = NEON_stab_plots_6yr.data$alpha_div)
head(alpha_spe_sta_6yr_newdat.lme)
alpha_spe_sta_6yr_newdat.lme$predlme = predict(alpha_spe_sta_6yr_lmm.fit, newdata = alpha_spe_sta_6yr_newdat.lme, level = 0)
alpha_spe_sta_6yr_des = model.matrix(formula(alpha_spe_sta_6yr_lmm.fit)[-2], alpha_spe_sta_6yr_newdat.lme)
alpha_spe_sta_6yr_predvar = diag( alpha_spe_sta_6yr_des %*% vcov(alpha_spe_sta_6yr_lmm.fit) %*% t(alpha_spe_sta_6yr_des))
alpha_spe_sta_6yr_newdat.lme$lower = with(alpha_spe_sta_6yr_newdat.lme, predlme - 2*sqrt(alpha_spe_sta_6yr_predvar))
alpha_spe_sta_6yr_newdat.lme$upper = with(alpha_spe_sta_6yr_newdat.lme, predlme + 2*sqrt(alpha_spe_sta_6yr_predvar))
################################################# plotting ######################################################
alpha_spe_sta_6yr.plot <- ggplot(NEON_stab_plots_6yr.data, aes(y = spe_sta, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  #geom_ribbon(data = alpha_spa_sta_6yr_newdat.lme, 
  #aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  #scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_spe_sta_6yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black", linetype ="dashed") +
  ylab(bquote(paste("Species stability ("*NULL[paste(log[10])]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "D") +
  annotate("text", x = 2.1, y = 1.0, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.008\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 2.1, y = 0.75, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.38\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

# 3.5 the diversity-stability relationships at population-level 
# fit the model examining the effects of alpha diversity on species asynchrony 
alpha_spe_asy_6yr_lmm.fit <- lme(spe_asy ~ alpha_div, random = ~1|siteID/nlcdClass, 
                                 correlation = corExp(form = ~ latitude + longitude|siteID/nlcdClass), data = NEON_stab_plots_6yr.data)
summary(alpha_spe_asy_6yr_lmm.fit)
anova(alpha_spe_asy_6yr_lmm.fit)
r.squaredGLMM(alpha_spe_asy_6yr_lmm.fit)
# test homoscedasticity
fligner.test(resid(alpha_spe_asy_6yr_lmm.fit), NEON_stab_plots_6yr.data$alpha_div)
# Assess normality of residuals
hist(alpha_spe_asy_6yr_lmm.fit$residuals)
# Inspect the model diagnostic metrics
plot(alpha_spe_asy_6yr_lmm.fit)

############################################## generating new data ###############################################
alpha_spe_asy_6yr_newdat.lme = data.frame(siteID = NEON_stab_plots_6yr.data$siteID,
                                          alpha_div = NEON_stab_plots_6yr.data$alpha_div)
head(alpha_spe_asy_6yr_newdat.lme)
alpha_spe_asy_6yr_newdat.lme$predlme = predict(alpha_spe_asy_6yr_lmm.fit, newdata = alpha_spe_asy_6yr_newdat.lme, level = 0)
alpha_spe_asy_6yr_des = model.matrix(formula(alpha_spe_asy_6yr_lmm.fit)[-2], alpha_spe_asy_6yr_newdat.lme)
alpha_spe_asy_6yr_predvar = diag( alpha_spe_asy_6yr_des %*% vcov(alpha_spe_asy_6yr_lmm.fit) %*% t(alpha_spe_asy_6yr_des))
alpha_spe_asy_6yr_newdat.lme$lower = with(alpha_spe_asy_6yr_newdat.lme, predlme - 2*sqrt(alpha_spe_asy_6yr_predvar))
alpha_spe_asy_6yr_newdat.lme$upper = with(alpha_spe_asy_6yr_newdat.lme, predlme + 2*sqrt(alpha_spe_asy_6yr_predvar))
################################################# plotting ######################################################
alpha_spe_asy_6yr.plot <- ggplot(NEON_stab_plots_6yr.data, aes(y = spe_asy, x = alpha_div, colour = factor(siteID))) +
  geom_point(size = 2.5, alpha = 0.3) +
  geom_smooth(aes(group = siteID), method = "lm", formula = y~x, se=F, size =0.5) +
  geom_ribbon(data = alpha_spe_asy_6yr_newdat.lme, 
              aes(y = NULL, ymin = lower, ymax = upper, color = NULL), alpha = .15) +  # add Confidence interval
  scale_fill_manual(values=c("gray")) + # fill the color
  geom_line(data = alpha_spe_asy_6yr_newdat.lme, aes(y = predlme), size = 1.5, color = "black") +
  ylab(bquote(paste("Species asynchrony ("*NULL[paste(log[10])]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "E") +
  annotate("text", x = 1.5, y = 2.1, size = 4, family = "Arial", 
           label ="paste(italic(R[m]^2), \" = 0.07 ***\")", parse=TRUE, hjust= 0) +
  annotate("text", x = 1.5, y = 1.8, size = 4, family = "Arial", 
           label ="paste(italic(R[c]^2), \" = 0.64 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  my_theme

