# Load packages
library(MuMIn)
library(ggplot2)

# load data
plot_dist.data <- read.csv('community_dis_dist.csv',header=T)
# check variables
variable.names(plot_dist.data)
##############################################################################################################
[1] "siteID" "plot_dist_km" "comdis" "plot_dist_km_ln" "comdis_ln" 
##############################################################################################################
##############################################################################################################
#################### transformation function
scaleFUN <- function(x) sprintf("%.1f", x)
############################################

site <- unique(plot_dist.data$siteID); site

all_dis.result <- c()
for (j in 1:length (site)) {
  data <- plot_dist.data[plot_dist.data$siteID==site[j], ]
  each.x <- data[ , 2] #spa_dist
  each.y <- data[ , 3] #community dissimilarity
  each.fit <- lm(log10(each.y) ~ log10(each.x))
    each.result <- c(unlist(summary(each.fit)$coefficients[1, 1:2]),
                     unlist(summary(each.fit)$coefficients[2, 1:4]), 
                     r2 <- r.squaredGLMM(each.fit))
    d <- c(as.character(site[j]), each.result)
    all_dis.result <- rbind(all_dis.result, d)
}
colnames(all_dis.result) <- c("siteID", "intercept", "se_intercept", "estimate", "se", "t value", "p", "r2m", "r2c")

## we have calculated the community dissimilarity between two pairwise plots 
## examine the relationships between community dissimilarity and the distance between two pairwise plots 

dist_all.plot <- ggplot(plot_dist.data, aes(y = log10(comdis), x = log10(plot_dist_km))) +
  geom_point(size = 1, alpha = 0.1, color = "gray30") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1) +
  facet_wrap(~ siteID, ncol = 6, scales = "free") +
  ylab(bquote(paste("Community dissimilarity ("*log[10]*"(Bray-Curtis index)"*") "))) +
  xlab(bquote(paste("Plot distance ("*log[10]*"(km)"*") "))) +
  #scale_y_continuous(labels=scaleFUN) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 16, family = "Arial", color = "black"),
        legend.text = element_text(size = 16, family = "Arial", color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_rect(colour = "black", size = 0.5),
        axis.text.y = element_text(size = 10, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 10, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(color = "black", size = 0.5),
        strip.text = element_text(size = 10, family = "Arial", color = "black"))

