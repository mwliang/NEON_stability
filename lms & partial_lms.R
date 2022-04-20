# Load packages
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
[17] "no.plot"                   "alpha_div"                 "beta_div1"                 "gamma_div"                
[21] "beta_div2"                 "tau_div"                   "spe_sta"                   "spe_asyn"                 
[25] "alpha_sta"                 "spa_asyn1"                 "gamma_sta"                 "spa_asyn2"                
[29] "tau_sta"                   "slope_alpha_sta.alpha_div" "slope_beta_beta"           "slope_gamma_gamma" 
##############################################################################################################
##############################################################################################################

############################################
# standardize function (z-score)
my_scale <- function(x){
  x1 <- (x-mean(x))/sd(x)
  return(x1)
}
############################################

############################################
#transformation function
scaleFUN <- function(x) sprintf("%.1f", x)
############################################

###############################################################################################
my_theme <- theme(legend.position = "none",
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
# 1 for >= 4-year observations (N=36)
NEON_stab_4yr.data <- NEON_stab.data
## 1.1 the results based on the simple linear regression

# 1.1.1 the diversity-stability relationships at alpha scale (Fig. 2A and Supplementary Table 1)
# fit the model examining the effects of alpha diversity on alpha stability
alpha_4yr.lm.fit <- lm(alpha_sta ~ alpha_div, data = NEON_stab_4yr.data)
summary(alpha_4yr.lm.fit)
# Assess normality of residuals
hist(alpha_4yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(alpha_4yr.lm.fit, which=1:4)

# compare with the nonlinear regression models (i.e., quadratic)
# alpha_4yr.qua.fit <- lm(alpha_sta ~ poly(alpha_div,2), data = NEON_stab_4yr.data)
# summary(alpha_4yr.qua.fit)
# anova(alpha_4yr.lm.fit, alpha_4yr.qua.fit)

alpha_4yr.plot <- ggplot(NEON_stab_4yr.data, aes(y = alpha_sta, x = alpha_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#0072B2", fill = "#0072B2", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#0072B2") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  theme_bw() +
  my_theme

# 1.1.2 the diversity-stability relationships at gamma scale (Fig. 2B and Supplementary Table 1)
# fit the model examining the effects of gamma diversity on gamma stability
gamma_4yr.lm.fit <- lm(gamma_sta ~ gamma_div, data = NEON_stab_4yr.data)
summary(gamma_4yr.lm.fit)
# Assess normality of residuals
hist(gamma_4yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(gamma_4yr.lm.fit, which=1:4)

# compare with the nonlinear regression models (i.e., quadratic)
#gamma_4yr.qua.fit <- lm(gamma_sta ~ poly(gamma_div,2), data = NEON_stab_4yr.data)
#summary(gamma_4yr.qua.fit)
#anova(gamma_4yr.lm.fit, gamma_4yr.qua.fit)

gamma_4yr.plot <- ggplot(NEON_stab_4yr.data, aes(y = gamma_sta, x = gamma_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#D55E00", fill = "#D55E00", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#D55E00") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  theme_bw() +
  my_theme

# 1.1.3 the diversity-stability relationships at tau scale (Fig. 2C and Supplementary Table 1)
# fit the model examining the effects of tau diversity on tau stability
tau_4yr.lm.fit <- lm(tau_sta ~ tau_div, data = NEON_stab_4yr.data)
summary(tau_4yr.lm.fit)
# Assess normality of residuals
hist(tau_4yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(tau_4yr.lm.fit, which=1:4)

# compare with the nonlinear regression models (i.e., quadratic)
#tau_4yr.qua.fit <- lm(tau_sta ~ poly(tau_div,2), data = NEON_stab_4yr.data)
#summary(tau_4yr.qua.fit)
#anova(tau_4yr.lm.fit, tau_4yr.qua.fit)

tau_4yr.plot <- ggplot(NEON_stab_4yr.data, aes(y = tau_sta, x = tau_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#009E73", fill = "#009E73", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#009E73") +
  ylab(bquote(paste(tau*" stability ("*NULL[paste(log[10]*"("*tau[S]*")")]*")"))) +
  xlab(bquote(paste(tau*" diversity ("*NULL[paste(log[10]*"("*tau[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  theme_bw() +
  my_theme

# 1.1.4 the diversity-stability relationships at beta1 scale (Fig. 2D and Supplementary Table 1)
# fit the model examining the effects of beta diversity 1 on spatial asynchrony 1 (beta stability 1)
beta1_4yr.lm.fit <- lm(spa_asyn1 ~ beta_div1, data = NEON_stab_4yr.data)
summary(beta1_4yr.lm.fit)
# Assess normality of residuals
hist(beta1_4yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta1_4yr.lm.fit, which=1:4)

# compare with the nonlinear regression models (i.e., quadratic)
#beta1_4yr.qua.fit <- lm(spa_asyn1 ~ poly(beta_div1,2), data = NEON_stab_4yr.data)
#summary(beta1_4yr.qua.fit)
#anova(beta1_4yr.lm.fit, beta1_4yr.qua.fit)

beta1_4yr.plot <- ggplot(NEON_stab_4yr.data, aes(y = spa_asyn1, x = beta_div1)) +
  geom_point(shape = 21, size = 3.5, colour = "gray30", fill = "gray30", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="gray30") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "D") +
  theme_bw() +
  my_theme

# 1.1.5 the diversity-stability relationships at beta2 scale (Fig. 2E and Supplementary Table 1)
# fit the model examining the effects of beta diversity 2 on spatial asynchrony 2 (beta stability 2)
beta2_4yr.lm.fit <- lm(spa_asyn2 ~ beta_div2, data = NEON_stab_4yr.data)
summary(beta2_4yr.lm.fit)
# Assess normality of residuals
hist(beta2_4yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta2_4yr.lm.fit, which=1:4)

# compare with the nonlinear regression models (i.e., quadratic)
#beta2_4yr.qua.fit <- lm(spa_asyn2 ~ poly(beta_div2,2), data = NEON_stab_4yr.data)
#summary(beta2_4yr.qua.fit)
#anova(beta2_4yr.lm.fit, beta2_4yr.qua.fit)

beta2_4yr.plot <- ggplot(NEON_stab_4yr.data, aes(y = spa_asyn2, x = beta_div2)) +
  geom_point(shape = 21, size = 3.5, colour = "black", fill = "black", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab(bquote(paste(beta^(gamma*"¡ú"*tau)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  xlab(bquote(paste(beta^(gamma*"¡ú"*tau)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "E") +
  theme_bw() +
  my_theme


## 1.2 the results based on the linear partial regression
# extract the residuals from the linear model,
# we used MAP and MAT because both MAP and MAT has been generally recognized 
# as the predominant environmental factors along the geographic gradients.
dd <- NEON_stab_4yr.data

values_4yr_plm.result <- c()
for (i in 19:30){
  each.depend <- dd[ , i]
  each.x1 <- dd[ , 10] # MAT_C
  each.x2 <- dd[ , 11] # MAP_mm
  each.fit <- lm(each.depend ~ my_scale(each.x1) + my_scale(each.x2))
  resdat.result <- each.fit$residuals
  values_plm.result <- cbind(values_4yr_plm.result, resdat.result)
}
colnames(values_4yr_plm.result) <- c("alpha_div_p",
                                     "beta_div1_p",                     
                                     "gamma_div_p",                
                                     "beta_div2_p",                     
                                     "tau_div_p",                   
                                     "spe_sta_p",                   
                                     "spe_asyn_p",                 
                                     "alpha_sta_p",                 
                                     "spa_asyn1_p",                 
                                     "gamma_sta_p",                 
                                     "spa_asyn2_p",                
                                     "tau_sta_p")
#write.csv(values_4yr_plm.result, file = "values_4yr_plm.result.csv")

NEON_stab_4yr_p.data <- as.data.frame(values_4yr_plm.result)

# 1.2.1 the diversity-stability relationships at alpha scale (Supplementary Fig. 1A and Table 2)
# fit the model examining the effects of alpha diversity on alpha stability
alpha_4yr_p.lm.fit <- lm(alpha_sta_p ~ alpha_div_p, data = NEON_stab_4yr_p.data)
summary(alpha_4yr_p.lm.fit)
# Assess normality of residuals
hist(alpha_4yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(alpha_4yr_p.lm.fit, which=1:4)

alpha_4yr_p.plot <- ggplot(NEON_stab_4yr_p.data, aes(y = alpha_sta_p, x = alpha_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#0072B2", fill = "#0072B2", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#0072B2") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  theme_bw() +
  my_theme

# 1.2.2 the diversity-stability relationships at gamma scale (Supplementary Fig. 1B and Table 2)
# fit the model examining the effects of gamma diversity on gamma stability
gamma_4yr_p.lm.fit <- lm(gamma_sta_p ~ gamma_div_p, data = NEON_stab_4yr_p.data)
summary(gamma_4yr_p.lm.fit)
# Assess normality of residuals
hist(gamma_4yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(gamma_4yr_p.lm.fit, which=1:4)

gamma_4yr_p.plot <- ggplot(NEON_stab_4yr_p.data, aes(y = gamma_sta_p, x = gamma_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#D55E00", fill = "#D55E00", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#D55E00") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  theme_bw() +
  my_theme

# 1.2.3 the diversity-stability relationships at tau scale (Supplementary Fig. 1C and Table 2)
# fit the model examining the effects of tau diversity on tau stability
tau_4yr_p.lm.fit <- lm(tau_sta_p ~ tau_div_p, data = NEON_stab_4yr_p.data)
summary(tau_4yr_p.lm.fit)
# Assess normality of residuals
hist(tau_4yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(tau_4yr_p.lm.fit, which=1:4)

tau_4yr_p.plot <- ggplot(NEON_stab_4yr_p.data, aes(y = tau_sta_p, x = tau_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#009E73", fill = "#009E73", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#009E73") +
  ylab(bquote(paste(tau*" stability ("*NULL[paste(log[10]*"("*tau[S]*")")]*")"))) +
  xlab(bquote(paste(tau*" diversity ("*NULL[paste(log[10]*"("*tau[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  theme_bw() +
  my_theme

# 1.2.4 the diversity-stability relationships at beta1 scale (Supplementary Fig. 2D and Table 1)
# fit the model examining the effects of of beta diversity 1 on spatial asynchrony (beta stability 1)
beta1_4yr_p.lm.fit <- lm(spa_asyn1_p ~ beta_div1_p, data = NEON_stab_4yr_p.data)
summary(beta1_4yr_p.lm.fit)
# Assess normality of residuals
hist(beta1_4yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta1_4yr_p.lm.fit, which=1:4)

beta1_4yr_p.plot <- ggplot(NEON_stab_4yr_p.data, aes(y = spa_asyn1_p, x = beta_div1_p)) +
  geom_point(shape = 21, size = 3.5, colour = "gray30", fill = "gray30", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="gray30") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "D") +
  theme_bw() +
  my_theme

# 1.2.5 the diversity-stability relationships at beta2 scale (Supplementary Fig. 1E and Table 2)
# fit the model examining the effects of of beta diversity 2 on spatial asynchrony (beta stability 2)
beta2_4yr_p.lm.fit <- lm(spa_asyn2_p ~ beta_div2_p, data = NEON_stab_4yr_p.data)
summary(beta2_4yr_p.lm.fit)
# Assess normality of residuals
hist(beta2_4yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta2_4yr_p.lm.fit, which=1:4)

beta2_4yr_p.plot <- ggplot(NEON_stab_4yr_p.data, aes(y = spa_asyn2_p, x = beta_div2_p)) +
  geom_point(shape = 21, size = 3.5, colour = "black", fill = "black", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab(bquote(paste(beta^(gamma*"¡ú"*tau)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  xlab(bquote(paste(beta^(gamma*"¡ú"*tau)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "E") +
  theme_bw() +
  my_theme



############################################
# 2 for >= 5-year observations (N=24)
NEON_stab_5yr.data <- subset(NEON_stab.data, NEON_stab.data$duration > 4)
## 2.1 the results based on the simple linear regression

# 2.1.1 the diversity-stability relationships at alpha scale (Extended Data Fig. 3A and Supplementary Table 1)
# fit the model examining the effects of alpha diversity on alpha stability
alpha_5yr.lm.fit <- lm(alpha_sta ~ alpha_div, data = NEON_stab_5yr.data)
summary(alpha_5yr.lm.fit)
# Assess normality of residuals
hist(alpha_5yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(alpha_5yr.lm.fit, which=1:4)

alpha_5yr.plot <- ggplot(NEON_stab_5yr.data, aes(y = alpha_sta, x = alpha_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#0072B2", fill = "#0072B2", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#0072B2") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  theme_bw() +
  my_theme

# 2.1.2 the diversity-stability relationships at gamma scale (Extended Data Fig. 3B and Supplementary Table 1)
# fit the model examining the effects of gamma diversity on gamma stability
gamma_5yr.lm.fit <- lm(gamma_sta ~ gamma_div, data = NEON_stab_5yr.data)
summary(gamma_5yr.lm.fit)
# Assess normality of residuals
hist(gamma_5yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(gamma_5yr.lm.fit, which=1:4)

gamma_5yr.plot <- ggplot(NEON_stab_5yr.data, aes(y = gamma_sta, x = gamma_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#D55E00", fill = "#D55E00", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#D55E00") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  theme_bw() +
  my_theme

# 2.1.3 the diversity-stability relationships at tau scale (Extended Data Fig. 3C and Supplementary Table 1)
# fit the model examining the effects of tau diversity on tau stability
tau_5yr.lm.fit <- lm(tau_sta ~ tau_div, data = NEON_stab_5yr.data)
summary(tau_5yr.lm.fit)
# Assess normality of residuals
hist(tau_5yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(tau_5yr.lm.fit, which=1:4)

tau_5yr.plot <- ggplot(NEON_stab_5yr.data, aes(y = tau_sta, x = tau_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#009E73", fill = "#009E73", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#009E73") +
  ylab(bquote(paste(tau*" stability ("*NULL[paste(log[10]*"("*tau[S]*")")]*")"))) +
  xlab(bquote(paste(tau*" diversity ("*NULL[paste(log[10]*"("*tau[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  theme_bw() +
  my_theme

# 2.1.4 the diversity-stability relationships at beta1 scale (Extended Data Fig. 3D and Supplementary Table 1)
# fit the model examining the effects of beta diversity 1 on spatial asynchrony 1 (beta stability 1)
beta1_5yr.lm.fit <- lm(spa_asyn1 ~ beta_div1, data = NEON_stab_5yr.data)
summary(beta1_5yr.lm.fit)
# Assess normality of residuals
hist(beta1_5yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta1_5yr.lm.fit, which=1:4)

beta1_5yr.plot <- ggplot(NEON_stab_5yr.data, aes(y = spa_asyn1, x = beta_div1)) +
  geom_point(shape = 21, size = 3.5, colour = "gray30", fill = "gray30", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="gray30") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "D") +
  theme_bw() +
  my_theme

# 2.1.5 the diversity-stability relationships at beta2 scale (Extended Data Fig. 3E and Supplementary Table 1)
# fit the model examining the effects of beta diversity 2 on spatial asynchrony 2 (beta stability 2)
beta2_5yr.lm.fit <- lm(spa_asyn2 ~ beta_div2, data = NEON_stab_5yr.data)
summary(beta2_5yr.lm.fit)
# Assess normality of residuals
hist(beta2_5yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta2_5yr.lm.fit, which=1:4)

beta2_5yr.plot <- ggplot(NEON_stab_5yr.data, aes(y = spa_asyn2, x = beta_div2)) +
  geom_point(shape = 21, size = 3.5, colour = "black", fill = "black", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab(bquote(paste(beta^(gamma*"¡ú"*tau)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  xlab(bquote(paste(beta^(gamma*"¡ú"*tau)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "E") +
  theme_bw() +
  my_theme


## 2.2 the results based on the linear partial regression
# extract the residuals from linear model first
# we used MAP and MAT because both MAP and MAT has been generally recognized 
# as the predominant environmental factors along the geographic gradients.
dd <- NEON_stab_5yr.data

values_5yr_plm.result <- c()
for (i in 19:30){
  each.depend <- dd[ , i]
  each.x1 <- dd[ , 10] # MAT_C
  each.x2 <- dd[ , 11] # MAP_mm
  each.fit <- lm(each.depend ~ my_scale(each.x1) + my_scale(each.x2))
  resdat.result <- each.fit$residuals
  values_plm.result <- cbind(values_5yr_plm.result, resdat.result)
}
colnames(values_5yr_plm.result) <- c("alpha_div_p",
                                     "beta_div1_p",                     
                                     "gamma_div_p",                
                                     "beta_div2_p",                     
                                     "tau_div_p",                   
                                     "spe_sta_p",                   
                                     "spe_asyn_p",                 
                                     "alpha_sta_p",                 
                                     "spa_asyn1_p",                 
                                     "gamma_sta_p",                 
                                     "spa_asyn2_p",                
                                     "tau_sta_p")
#write.csv(values_5yr_plm.result, file = "values_5yr_plm.result.csv")

NEON_stab_5yr_p.data <- as.data.frame(values_5yr_plm.result)

# 2.2.1 the diversity-stability relationships at alpha scale (Supplementary Table 2)
# fit the model examining the effects of alpha diversity on alpha stability
alpha_5yr_p.lm.fit <- lm(alpha_sta_p ~ alpha_div_p, data = NEON_stab_5yr_p.data)
summary(alpha_5yr_p.lm.fit)
# Assess normality of residuals
hist(alpha_5yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(alpha_5yr_p.lm.fit, which=1:4)

alpha_5yr_p.plot <- ggplot(NEON_stab_5yr_p.data, aes(y = alpha_sta_p, x = alpha_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#0072B2", fill = "#0072B2", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#0072B2") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  theme_bw() +
  my_theme

# 2.2.2 the diversity-stability relationships at gamma scale (Supplementary Table 2)
# fit the model examining the effects of gamma diversity on gamma stability
gamma_5yr_p.lm.fit <- lm(gamma_sta_p ~ gamma_div_p, data = NEON_stab_5yr_p.data)
summary(gamma_5yr_p.lm.fit)
# Assess normality of residuals
hist(gamma_5yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(gamma_5yr_p.lm.fit, which=1:4)

gamma_5yr_p.plot <- ggplot(NEON_stab_5yr_p.data, aes(y = gamma_sta_p, x = gamma_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#D55E00", fill = "#D55E00", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#D55E00") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  theme_bw() +
  my_theme

# 2.2.3 the diversity-stability relationships at tau scale (Supplementary Table 2)
# fit the model examining the effects of tau diversity on tau stability
tau_5yr_p.lm.fit <- lm(tau_sta_p ~ tau_div_p, data = NEON_stab_5yr_p.data)
summary(tau_5yr_p.lm.fit)
# Assess normality of residuals
hist(tau_5yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(tau_5yr_p.lm.fit, which=1:4)

tau_5yr_p.plot <- ggplot(NEON_stab_5yr_p.data, aes(y = tau_sta_p, x = tau_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#009E73", fill = "#009E73", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#009E73") +
  ylab(bquote(paste(tau*" stability ("*NULL[paste(log[10]*"("*tau[S]*")")]*")"))) +
  xlab(bquote(paste(tau*" diversity ("*NULL[paste(log[10]*"("*tau[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  theme_bw() +
  my_theme

# 2.2.4 the diversity-stability relationships at beta1 scale (Supplementary Table 2)
# fit the model examining the effects of of beta diversity 1 on spatial asynchrony (beta stability 1)
beta1_5yr_p.lm.fit <- lm(spa_asyn1_p ~ beta_div1_p, data = NEON_stab_5yr_p.data)
summary(beta1_5yr_p.lm.fit)
# Assess normality of residuals
hist(beta1_5yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta1_5yr_p.lm.fit, which=1:4)

beta1_5yr_p.plot <- ggplot(NEON_stab_5yr_p.data, aes(y = spa_asyn1_p, x = beta_div1_p)) +
  geom_point(shape = 21, size = 3.5, colour = "gray30", fill = "gray30", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="gray30") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "D") +
  theme_bw() +
  my_theme

# 2.2.5 the diversity-stability relationships at beta2 scale (Supplementary Table 2)
# fit the model examining the effects of of beta diversity 2 on spatial asynchrony (beta stability 2)
beta2_5yr_p.lm.fit <- lm(spa_asyn2_p ~ beta_div2_p, data = NEON_stab_5yr_p.data)
summary(beta2_5yr_p.lm.fit)
# Assess normality of residuals
hist(beta2_5yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta2_5yr_p.lm.fit, which=1:4)

beta2_5yr_p.plot <- ggplot(NEON_stab_5yr_p.data, aes(y = spa_asyn2_p, x = beta_div2_p)) +
  geom_point(shape = 21, size = 3.5, colour = "black", fill = "black", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab(bquote(paste(beta^(gamma*"¡ú"*tau)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  xlab(bquote(paste(beta^(gamma*"¡ú"*tau)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "E") +
  theme_bw() +
  my_theme


############################################
# 3 for >= 6-year observations (N=14)
NEON_stab_6yr.data <- subset(NEON_stab.data, NEON_stab.data$duration > 5)
## 3.1 the results based on the simple linear regression

# 3.1.1 the diversity-stability relationships at alpha scale (Extended Data Fig. 4A and Supplementary Table 1)
# fit the model examining the effects of alpha diversity on alpha stability
alpha_6yr.lm.fit <- lm(alpha_sta ~ alpha_div, data = NEON_stab_6yr.data)
summary(alpha_6yr.lm.fit)
# Assess normality of residuals
hist(alpha_6yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(alpha_6yr.lm.fit, which=1:4)

alpha_6yr.plot <- ggplot(NEON_stab_6yr.data, aes(y = alpha_sta, x = alpha_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#0072B2", fill = "#0072B2", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#0072B2") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  theme_bw() +
  my_theme

# 3.1.2 the diversity-stability relationships at gamma scale (Extended Data Fig. 4B and Supplementary Table 1)
# fit the model examining the effects of gamma diversity on gamma stability
gamma_6yr.lm.fit <- lm(gamma_sta ~ gamma_div, data = NEON_stab_6yr.data)
summary(gamma_6yr.lm.fit)
# Assess normality of residuals
hist(gamma_6yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(gamma_6yr.lm.fit, which=1:4)

gamma_6yr.plot <- ggplot(NEON_stab_6yr.data, aes(y = gamma_sta, x = gamma_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#D55E00", fill = "#D55E00", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#D55E00") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  theme_bw() +
  my_theme

# 3.1.3 the diversity-stability relationships at tau scale (Extended Data Fig. 4C and Supplementary Table 1)
# fit the model examining the effects of tau diversity on tau stability
tau_6yr.lm.fit <- lm(tau_sta ~ tau_div, data = NEON_stab_6yr.data)
summary(tau_6yr.lm.fit)
# Assess normality of residuals
hist(tau_6yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(tau_6yr.lm.fit, which=1:4)

tau_6yr.plot <- ggplot(NEON_stab_6yr.data, aes(y = tau_sta, x = tau_div)) +
  geom_point(shape = 21, size = 3.5, colour = "#009E73", fill = "#009E73", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#009E73") +
  ylab(bquote(paste(tau*" stability ("*NULL[paste(log[10]*"("*tau[S]*")")]*")"))) +
  xlab(bquote(paste(tau*" diversity ("*NULL[paste(log[10]*"("*tau[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  theme_bw() +
  my_theme

# 3.1.4 the diversity-stability relationships at beta1 scale (Extended Data Fig. 4D and Supplementary Table 1)
# fit the model examining the effects of beta diversity 1 on spatial asynchrony 1 (beta stability 1)
beta1_6yr.lm.fit <- lm(spa_asyn1 ~ beta_div1, data = NEON_stab_6yr.data)
summary(beta1_6yr.lm.fit)
# Assess normality of residuals
hist(beta1_6yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta1_6yr.lm.fit, which=1:4)

beta1_6yr.plot <- ggplot(NEON_stab_6yr.data, aes(y = spa_asyn1, x = beta_div1)) +
  geom_point(shape = 21, size = 3.5, colour = "gray30", fill = "gray30", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="gray30") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "D") +
  theme_bw() +
  my_theme

# 3.1.5 the diversity-stability relationships at beta2 scale (Extended Data Fig. 4E and Supplementary Table 1)
# fit the model examining the effects of beta diversity 2 on spatial asynchrony 2 (beta stability 2)
beta2_6yr.lm.fit <- lm(spa_asyn2 ~ beta_div2, data = NEON_stab_6yr.data)
summary(beta2_6yr.lm.fit)
# Assess normality of residuals
hist(beta2_6yr.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta2_6yr.lm.fit, which=1:4)

beta2_6yr.plot <- ggplot(NEON_stab_6yr.data, aes(y = spa_asyn2, x = beta_div2)) +
  geom_point(shape = 21, size = 3.5, colour = "black", fill = "black", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black", linetype ="dashed", se=F) +
  ylab(bquote(paste(beta^(gamma*"¡ú"*tau)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  xlab(bquote(paste(beta^(gamma*"¡ú"*tau)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "E") +
  theme_bw() +
  my_theme


## 3.2 the results based on the linear partial regression
# extract the residuals from linear model first
# we used MAP and MAT because both MAP and MAT has been generally recognized 
# as the predominant environmental factors along the geographic gradients.
dd <- NEON_stab_6yr.data

values_6yr_plm.result <- c()
for (i in 19:30){
  each.depend <- dd[ , i]
  each.x1 <- dd[ , 10] # MAT_C
  each.x2 <- dd[ , 11] # MAP_mm
  each.fit <- lm(each.depend ~ my_scale(each.x1) + my_scale(each.x2))
  resdat.result <- each.fit$residuals
  values_plm.result <- cbind(values_6yr_plm.result, resdat.result)
}
colnames(values_6yr_plm.result) <- c("alpha_div_p",
                                     "beta_div1_p",                     
                                     "gamma_div_p",                
                                     "beta_div2_p",                     
                                     "tau_div_p",                   
                                     "spe_sta_p",                   
                                     "spe_asyn_p",                 
                                     "alpha_sta_p",                 
                                     "spa_asyn1_p",                 
                                     "gamma_sta_p",                 
                                     "spa_asyn2_p",                
                                     "tau_sta_p")
#write.csv(values_6yr_plm.result, file = "values_6yr_plm.result.csv")

NEON_stab_6yr_p.data <- as.data.frame(values_6yr_plm.result)

# 3.2.1 the diversity-stability relationships at alpha scale (Supplementary Table 2)
# fit the model examining the effects of alpha diversity on alpha stability
alpha_6yr_p.lm.fit <- lm(alpha_sta_p ~ alpha_div_p, data = NEON_stab_6yr_p.data)
summary(alpha_6yr_p.lm.fit)
# Assess normality of residuals
hist(alpha_6yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(alpha_6yr_p.lm.fit, which=1:4)

alpha_6yr_p.plot <- ggplot(NEON_stab_6yr_p.data, aes(y = alpha_sta_p, x = alpha_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#0072B2", fill = "#0072B2", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#0072B2") +
  ylab(bquote(paste(alpha*" stability ("*NULL[paste(log[10]*"("*alpha[S]*")")]*")"))) +
  xlab(bquote(paste(alpha*" diversity ("*NULL[paste(log[10]*"("*alpha[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "A") +
  theme_bw() +
  my_theme

# 3.2.2 the diversity-stability relationships at gamma scale (Supplementary Table 2)
# fit the model examining the effects of gamma diversity on gamma stability
gamma_6yr_p.lm.fit <- lm(gamma_sta_p ~ gamma_div_p, data = NEON_stab_6yr_p.data)
summary(gamma_6yr_p.lm.fit)
# Assess normality of residuals
hist(gamma_6yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(gamma_6yr_p.lm.fit, which=1:4)

gamma_6yr_p.plot <- ggplot(NEON_stab_6yr_p.data, aes(y = gamma_sta_p, x = gamma_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#D55E00", fill = "#D55E00", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#D55E00") +
  ylab(bquote(paste(gamma*" stability ("*NULL[paste(log[10]*"("*gamma[S]*")")]*")"))) +
  xlab(bquote(paste(gamma*" diversity ("*NULL[paste(log[10]*"("*gamma[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "B") +
  theme_bw() +
  my_theme

# 3.2.3 the diversity-stability relationships at tau scale (Supplementary Table 2)
# fit the model examining the effects of tau diversity on tau stability
tau_6yr_p.lm.fit <- lm(tau_sta_p ~ tau_div_p, data = NEON_stab_6yr_p.data)
summary(tau_6yr_p.lm.fit)
# Assess normality of residuals
hist(tau_6yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(tau_6yr_p.lm.fit, which=1:4)

tau_6yr_p.plot <- ggplot(NEON_stab_6yr_p.data, aes(y = tau_sta_p, x = tau_div_p)) +
  geom_point(shape = 21, size = 3.5, colour = "#009E73", fill = "#009E73", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="#009E73") +
  ylab(bquote(paste(tau*" stability ("*NULL[paste(log[10]*"("*tau[S]*")")]*")"))) +
  xlab(bquote(paste(tau*" diversity ("*NULL[paste(log[10]*"("*tau[D]*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "C") +
  theme_bw() +
  my_theme

# 3.2.4 the diversity-stability relationships at beta1 scale (Supplementary Table 2)
# fit the model examining the effects of of beta diversity 1 on spatial asynchrony (beta stability 1)
beta1_6yr_p.lm.fit <- lm(spa_asyn1_p ~ beta_div1_p, data = NEON_stab_6yr_p.data)
summary(beta1_6yr_p.lm.fit)
# Assess normality of residuals
hist(beta1_6yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta1_6yr_p.lm.fit, which=1:4)

beta1_6yr_p.plot <- ggplot(NEON_stab_6yr_p.data, aes(y = spa_asyn1_p, x = beta_div1_p)) +
  geom_point(shape = 21, size = 3.5, colour = "gray30", fill = "gray30", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="gray30") +
  ylab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  xlab(bquote(paste(beta^(alpha*"¡ú"*gamma)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(alpha*"¡ú"*gamma)*")")]*")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "D") +
  theme_bw() +
  my_theme

# 3.2.5 the diversity-stability relationships at beta2 scale (Supplementary Table 2)
# fit the model examining the effects of of beta diversity 2 on spatial asynchrony (beta stability 2)
beta2_6yr_p.lm.fit <- lm(spa_asyn2_p ~ beta_div2_p, data = NEON_stab_6yr_p.data)
summary(beta2_6yr_p.lm.fit)
# Assess normality of residuals
hist(beta2_6yr_p.lm.fit$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(beta2_6yr_p.lm.fit, which=1:4)

beta2_6yr_p.plot <- ggplot(NEON_stab_6yr_p.data, aes(y = spa_asyn2_p, x = beta_div2_p)) +
  geom_point(shape = 21, size = 3.5, colour = "black", fill = "black", alpha = 0.5) +
  #geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), size = 1, colour ="red", se=F) +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab(bquote(paste(beta^(gamma*"¡ú"*tau)*" stability ("*NULL[paste(log[10]*"("*beta[~S]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  xlab(bquote(paste(beta^(gamma*"¡ú"*tau)*" diversity ("*NULL[paste(log[10]*"("*beta[~D]^(gamma*"¡ú"*tau)*")")] * ")"))) +
  scale_y_continuous(labels=scaleFUN) +
  labs(tag = "E") +
  theme_bw() +
  my_theme