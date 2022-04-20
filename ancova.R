# Load packages
library(ggplot2)
library(ggiraphExtra)
library(ggExtra)
library(ggpubr)

# load data
ANCOVA.data <- read.csv('ANCOVA_div_sta.csv',header=T)
# check variables
variable.names(ANCOVA.data)
##############################################################################################################
[1] "state"        "biome_types"  "latitude"     "longitude"    "elevation_m" 
[6] "MAT_C"        "MAP_mm"       "Aridity"      "Temp_s"       "Prec_s"      
[11] "started_year" "end_year"     "duration"     "N"            "scale"       
[16] "scale1"       "div"          "sta"          "div_p"        "sta_p"
##############################################################################################################
############################################
# 1 for >= 4-year observations (N=36)
ANCOVA_4yr.data <- ANCOVA.data
## 1.1 the results based on the simple linear regression
# plot all scatter plots together with model information

ggscatter(ANCOVA_4yr.data, x = "div", y = "sta", color = "scale", add = "reg.line") +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = scale)) +
  scale_color_manual(values=c('#0072B2','#E69F00','#009E73','gray40',"black"),
                     breaks = c("alpha","gamma","tau","beta1","beta2")) +
  ylab(bquote(paste("Stability metrics ( "*log[10]*")"))) +
  xlab(bquote(paste("Diversity metrics ( "*log[10]*")"))) +
  labs(color = "Scale")

# "a" denotes alpha, gamma, and tau scale
# "b" denotes beta1 and beta2 scale

ANCOVA_4yr.data$scale1 <- factor(ANCOVA_4yr.data$scale1,levels=c("a","b"))

# Fig 2F left panel
# 1.1.1 comparing the effects of scale on the diversity-stability relationships at alpha, gamma, and tau scales
# select the data only for the alpha, gamma, and tau scales
ANCOVA_rich_4yr.data <- ANCOVA_4yr.data[!ANCOVA_4yr.data$scale1%in%c("b"),]
# fit the model examining the effects of scale
com_rich_4yr.aov <- aov(sta ~ div*scale, data = ANCOVA_rich_4yr.data)
summary(com_rich_4yr.aov)
# multiple comparisons between alpha vs. gamma, alpha vs. tau, and gamma vs. tau
TukeyHSD(com_rich_4yr.aov,'scale') 
# Assess normality of residuals
hist(com_rich_4yr.aov$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(com_rich_4yr.aov, which=1:4)

# Fig 2F right panel
# 1.1.2 comparing the effects of scale on the beta_div-spatial_asy relationships at beta1 and beta2 scales
# select the data only for the beta1 and beta2 scales
ANCOVA_beta_4yr.data <- ANCOVA_4yr.data[!ANCOVA_4yr.data$scale1%in%c("a"),]
# fit the model examining the effects of scale
asy_beta_4yr.aov <- aov(sta ~ div*scale, data = ANCOVA_beta_4yr.data)
summary(asy_beta_4yr.aov)
# Assess normality of residuals
hist(asy_beta_4yr.aov$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(asy_beta_4yr.aov, which=1:4)
dev.off()

## 1.2 the results based on the partial linear regression
# div_p and sta_p are the residuals of from the linear model taking account for
# the predominant environmental factors (MAP and MAT ) 

# plot all scatter plots together with model information
ggscatter(ANCOVA_4yr.data, x = "div_p", y = "sta_p", color = "scale", add = "reg.line") +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = scale)) +
  scale_color_manual(values=c('#0072B2','#E69F00','#009E73','gray40',"black"),
                     breaks = c("alpha","gamma","tau","beta1","beta2")) +
  ylab(bquote(paste("Stability metrics ( "*log[10]*")"))) +
  xlab(bquote(paste("Diversity metrics ( "*log[10]*")"))) +
  labs(color = "Scale")

# "a" denotes alpha, gamma, and tau scale
# "b" denotes beta1 and beta2 scale

ANCOVA_4yr.data$scale1 <- factor(ANCOVA_4yr.data$scale1,levels=c("a","b"))

# Supplementary Fig 2F left panel
# 1.2.1 comparing the effects of scale on the diversity-stability relationships at alpha, gamma, and tau scales
# select the data only for the alpha, gamma, and tau scales
ANCOVA_rich_p_4yr.data <- ANCOVA_4yr.data[!ANCOVA_4yr.data$scale1%in%c("b"),]
# fit the model examining the effects of scale
com_rich_p_4yr.aov <- aov(sta_p ~ div_p*scale, data = ANCOVA_rich_p_4yr.data)
summary(com_rich_p_4yr.aov)
# multiple comparisons between alpha vs. gamma, alpha vs. tau, and gamma vs. tau
TukeyHSD(com_rich_p_4yr.aov,'scale') 
# Assess normality of residuals
hist(com_rich_p_4yr.aov$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(com_rich_p_4yr.aov, which=1:4)

# Supplementary Fig 2F right panel
# 1.2.2 comparing the effects of scale on the beta_div-spatial_asy relationships at beta1 and beta2 scales
# select the data only for the beta1 and beta2 scales
ANCOVA_beta_p_4yr.data <- ANCOVA_4yr.data[!ANCOVA_4yr.data$scale1%in%c("a"),]
# fit the model examining the effects of scale
asy_beta_p_4yr.aov <- aov(sta_p ~ div_p*scale, data = ANCOVA_beta_p_4yr.data)
summary(asy_beta_p_4yr.aov)
# Assess normality of residuals
hist(asy_beta_p_4yr.aov$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(asy_beta_p_4yr.aov, which=1:4)



############################################
# 2 for >= 5-year observations (N=24)
ANCOVA_5yr.data <- subset(ANCOVA.data, ANCOVA.data$duration > 4)

## 2.1 the results based on the simple linear regression
# plot all scatter plots together with model information
ggscatter(ANCOVA_5yr.data, x = "div", y = "sta", color = "scale", add = "reg.line") +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = scale)) +
  scale_color_manual(values=c('#0072B2','#E69F00','#009E73','gray40',"black"),
                     breaks = c("alpha","gamma","tau","beta1","beta2")) +
  ylab(bquote(paste("Stability metrics ( "*log[10]*")"))) +
  xlab(bquote(paste("Diversity metrics ( "*log[10]*")"))) +
  labs(color = "Scale")

# "a" denotes alpha, gamma, and tau scale
# "b" denotes beta1 and beta2 scale

ANCOVA_5yr.data$scale1 <- factor(ANCOVA_5yr.data$scale1,levels=c("a","b"))

# Extended Data Figs. 7F left panel
# 2.1.1 comparing the effects of scale on the diversity-stability relationships at alpha, gamma, and tau scales
# select the data only for the alpha, gamma, and tau scales
ANCOVA_rich_5yr.data <- ANCOVA_5yr.data[!ANCOVA_5yr.data$scale1%in%c("b"),]
# fit the model examining the effects of scale
com_rich_5yr.aov <- aov(sta ~ div*scale, data = ANCOVA_rich_5yr.data)
summary(com_rich_5yr.aov)
# multiple comparisons between alpha vs. gamma, alpha vs. tau, and gamma vs. tau
TukeyHSD(com_rich_5yr.aov,'scale') 
# Assess normality of residuals
hist(com_rich_5yr.aov$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(com_rich_5yr.aov, which=1:4)

# Extended Data Figs. 7F right panel
# 2.1.2 comparing the effects of scale on the beta_div-spatial_asy relationships at beta1 and beta2 scales
# select the data only for the beta1 and beta2 scales
ANCOVA_beta_5yr.data <- ANCOVA_5yr.data[!ANCOVA_5yr.data$scale1%in%c("a"),]
# fit the model examining the effects of scale
asy_beta_5yr.aov <- aov(sta ~ div*scale, data = ANCOVA_beta_5yr.data)
summary(asy_beta_5yr.aov)
# Assess normality of residuals
hist(asy_beta_5yr.aov$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(asy_beta_5yr.aov, which=1:4)
dev.off()


############################################
# 3 for >= 5-year observations (N=14)
ANCOVA_6yr.data <- subset(ANCOVA.data, ANCOVA.data$duration > 5)

## 3.1 the results based on the simple linear regression
# plot all scatter plots together with model information
ggscatter(ANCOVA_6yr.data, x = "div", y = "sta", color = "scale", add = "reg.line") +
  stat_regline_equation(aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = scale)) +
  scale_color_manual(values=c('#0072B2','#E69F00','#009E73','gray40',"black"),
                     breaks = c("alpha","gamma","tau","beta1","beta2")) +
  ylab(bquote(paste("Stability metrics ( "*log[10]*")"))) +
  xlab(bquote(paste("Diversity metrics ( "*log[10]*")"))) +
  labs(color = "Scale")

# "a" denotes alpha, gamma, and tau scale
# "b" denotes beta1 and beta2 scale

ANCOVA_6yr.data$scale1 <- factor(ANCOVA_6yr.data$scale1,levels=c("a","b"))

# Extended Data Figs. 8F left panel
# 3.1.1 comparing the effects of scale on the diversity-stability relationships at alpha, gamma, and tau scales
# select the data only for the alpha, gamma, and tau scales
ANCOVA_rich_6yr.data <- ANCOVA_6yr.data[!ANCOVA_6yr.data$scale1%in%c("b"),]
# fit the model examining the effects of scale
com_rich_6yr.aov <- aov(sta ~ div*scale, data = ANCOVA_rich_6yr.data)
summary(com_rich_6yr.aov)
# multiple comparisons between alpha vs. gamma, alpha vs. tau, and gamma vs. tau
TukeyHSD(com_rich_6yr.aov,'scale') 
# Assess normality of residuals
hist(com_rich_6yr.aov$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(com_rich_6yr.aov, which=1:4)

# Extended Data Figs. 8F right panel
# 3.1.2 comparing the effects of scale on the beta_div-spatial_asy relationships at beta1 and beta2 scales
# select the data only for the beta1 and beta2 scales
ANCOVA_beta_6yr.data <- ANCOVA_6yr.data[!ANCOVA_6yr.data$scale1%in%c("a"),]
# fit the model examining the effects of scale
asy_beta_6yr.aov <- aov(sta ~ div*scale, data = ANCOVA_beta_6yr.data)
summary(asy_beta_6yr.aov)
# Assess normality of residuals
hist(asy_beta_6yr.aov$residuals)
# Inspect the model diagnostic metrics
par(mfrow=c(2,2))
plot(asy_beta_6yr.aov, which=1:4)
dev.off()
