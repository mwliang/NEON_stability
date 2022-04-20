
# Load packages
library(ggplot2)
library(cowplot)
library(sf)
library(ggspatial)
library(rnaturalearth)
library(rnaturalearthdata)
library(tools)
library(maps)
library(ggrepel)

# load data
NEON_stab.data <- read.csv('NEON_stab.csv',header=T)
# check variables
variable.names(NEON_stab.data)
######################################################################################################
[1] "no."                       "domain"                    "domain_name"               "siteID"                   
[5] "site_name"                 "state"                     "latitude"                  "longitude"                
[9] "elevation_m"               "MAT_C"                     "MAP_mm"                    "Temp_s"                   
[13] "Prec_s"                    "started_year"              "end_year"                  "duration"                 
[17] "no.plot"                   "spa_dist"                  "alpha_div"                 "beta_div1"                
[21] "gamma_div"                 "beta_div2"                 "tau_div"                   "spe_sta"                  
[25] "spe_asyn"                  "alpha_sta"                 "spa_asyn1"                 "gamma_sta"                
[29] "spa_asyn2"                 "tau_sta"                   "slope_alpha_sta.alpha_div" "slope_spa_asyn1.beta_div1"
[33] "slope_gamma_sta.gamma_div" "se_alpha"                  "se_beta1"                  "se_gamma"   
######################################################################################################

# Extended Data Fig. 1A 
states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
head(states)
states <- cbind(states, st_coordinates(st_centroid(states)))
states$ID <- toTitleCase(states$ID)
head(states)

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
usa <- subset(world, admin == "United States of America")

NEON_sites.matp <- ggplot() +
  geom_sf(data = world) +
  geom_sf(data = usa, fill = "antiquewhite1") +
  geom_sf(data = states, fill = NA, size = 0.2) +
  coord_sf(xlim = c(-170, -50), ylim = c(16, 75), expand = FALSE) +
  annotation_scale(location = "bl", 
                   pad_x = unit(0.3, "in"), 
                   pad_y = unit(0.9, "in"),
                   height = unit(0.1, "cm"), 
                   width_hint = 0.2) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.3, "in"), 
                         pad_y = unit(1.0, "in"),
                         style = north_arrow_fancy_orienteering) +
  geom_point(data = NEON_stab.data, aes(x = longitude, y = latitude, color = MAP_mm/10, size = MAT_C), alpha = 0.6) +
  scale_color_gradient(low = "orange", high = "blue") +
  guides(size=guide_legend(bquote(paste(MAT*" ("^o*C*")")))) +
  geom_text_repel(data = NEON_stab.data, aes(x = longitude, y = latitude, label = siteID),
                  size = 3, color = "grey30",
                  box.padding = unit(0.2, "lines"),
                  point.padding = unit(1, "lines"),
                  segment.color = 'grey50',
                  segment.size = 0.2,
                  # Draw an arrow from the label to the data point.
                  arrow = arrow(length = unit(0.01, 'npc')),
                  # Strength of the repulsion force.
                  force = 1,
                  # Maximum iterations of the naive repulsion algorithm O(n^2).
                  max.iter = 3e3) +
  labs(tag = "A", color = "MAP (cm)") +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme(#legend.position = c (0.15,0.5),
    #legend.position = "bottom",
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.text = element_text(size = 10, family = "Arial", color = "black"),
    panel.grid.major = element_line(colour = gray(0.5), linetype = "dashed", size = 0.2), 
    panel.background = element_rect(fill = "aliceblue"), 
    panel.border = element_rect(fill = NA),
    axis.text.y = element_text(size = 12, family = "Arial", color = "black"),
    axis.text.x = element_text(size = 12, family = "Arial", color = "black"),
    axis.title.y = element_text(size = 18, family = "Arial", color = "black"),
    axis.title.x = element_text(size = 18, family = "Arial", color = "black"),
    plot.tag = element_text(size = 22, family = "Arial", color = "black"))

########################################################################################
# Extended Data Fig. 1B
#devtools::install_github("valentinitnelav/plotbiomes")
library(plotbiomes)
#whittaker_base_plot(color_palette = NULL)

Whittaker_biomes$biome <- factor(Whittaker_biomes$biome, 
                                 levels = c("Tundra","Boreal forest", "Temperate seasonal forest",
                                            "Temperate rain forest", "Tropical rain forest", 
                                            "Tropical seasonal forest/savanna", "Subtropical desert",
                                            "Temperate grassland/desert", "Woodland/shrubland"))

mapt_biomes.plot <- ggplot() +
  geom_polygon(data = Whittaker_biomes, aes(x = temp_c, y = precp_cm, fill = biome),
               colour = "gray98", # colour of polygon border
               size = 0.5,
               alpha = 0.5) +    # thickness of polygon border
  # fill polygons with predefined colors (as in Ricklefs, 2008)
  scale_fill_manual(name   = "Whittaker biomes",
                    breaks = names(Ricklefs_colors),
                    labels = names(Ricklefs_colors),
                    values = Ricklefs_colors) +
  geom_point(data = NEON_stab.data, aes(x = MAT_C, y = MAP_mm/10), size = 3, color = "black", alpha = 0.5) +
  geom_text_repel(data = NEON_stab.data, aes(x = MAT_C, y = MAP_mm/10, label = siteID),
                  size = 3, color = "grey30",
                  box.padding = unit(0.2, "lines"),
                  point.padding = unit(1, "lines"),
                  segment.color = 'grey50',
                  segment.size = 0.2,
                  # Draw an arrow from the label to the data point.
                  arrow = arrow(length = unit(0.01, 'npc')),
                  # Strength of the repulsion force.
                  force = 1,
                  # Maximum iterations of the naive repulsion algorithm O(n^2).
                  max.iter = 3e3) +
  xlab(bquote(paste(MAT*" ("^o*C*")"))) +
  ylab("MAP (cm)") +
  labs(tag = "B") +
  #guides(size=guide_legend("Plant richness (site-level)")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.text = element_text(size = 16, family = "Arial", color = "black"),
        axis.title = element_text(size = 18, family = "Arial", color = "black"),
        plot.tag = element_text(size = 22, family = "Arial", color = "black"))






#NEON_stab.data$MAT_C <- as.numeric(NEON_stab.data$MAT_C)
#NEON_stab.data$MAP_mm <- as.numeric(NEON_stab.data$MAP_mm)

# Extended Data Fig. 1C
summary(lm(MAT_C ~ latitude, data = NEON_stab.data))

mat_lat.plot <- ggplot(NEON_stab.data, aes(y = MAT_C, x = latitude, colour = MAP_mm/10)) +
  geom_point(aes(size = MAT_C), alpha = 0.5) +
  scale_color_gradient(low = "orange", high = "blue") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab(bquote(paste(MAT*" ("^o*C*")"))) +
  xlab("Latitude") +
  labs(tag = "C") +
  ylim(-20, 30) +
  annotate("text", x = 40, y = 25, size = 5, family = "Arial", 
           label ="paste(italic(R^2), \" = 0.84 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(color = "black", size = 0.5),
        plot.tag = element_text(size = 20, family = "Arial", color = "black"))


# Extended Data Fig. 1D
summary(lm(MAP_mm/10 ~ latitude, data = NEON_stab.data))

map_lat.plot <- ggplot(NEON_stab.data, aes(y = MAP_mm/10, x = latitude, colour = MAP_mm/10)) +
  geom_point(aes(size = MAT_C), alpha = 0.5) +
  scale_color_gradient(low = "orange", high = "blue") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab("MAP (cm)") +
  xlab("Latitude") +
  labs(tag = "D") +
  annotate("text", x = 20, y = 230, size = 5, family = "Arial", 
           label ="paste(italic(R^2), \" = 0.12 *\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(color = "black", size = 0.5),
        plot.tag = element_text(size = 20, family = "Arial", color = "black"))


# Extended Data Fig. 1E
summary(lm(MAP_mm/10 ~ MAT_C, data = NEON_stab.data))

map_mat.plot <- ggplot(NEON_stab.data, aes(y = MAP_mm/10, x = latitude, colour = MAP_mm/10)) +
  geom_point(aes(size = MAT_C), alpha = 0.5) +
  scale_color_gradient(low = "orange", high = "blue") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab("MAP (cm)") +
  xlab(bquote(paste(MAT*" ("^o*C*")"))) +
  labs(tag = "E") +
  annotate("text", x = 20, y = 230, size = 5, family = "Arial", 
           label ="paste(italic(R^2), \" = 0.16 *\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(color = "black", size = 0.5),
        plot.tag = element_text(size = 20, family = "Arial", color = "black"))


########################################################################################
# Supplementary Fig. 5A
summary(lm(Temp_s ~ MAT_C, data = NEON_stab.data))

Temp_s_mat.plot <- ggplot(NEON_stab.data, aes(y = Temp_s, x = MAT_C, colour = MAP_mm)) +
  geom_point(aes(size = MAT_C), alpha = 0.5) +
  scale_color_gradient(low = "orange", high = "blue") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab(bquote(atop(paste("Temperature seasonality"), (paste("sd*100"))))) +
  xlab(bquote(paste(MAT*" ("^o*C*")"))) +
  labs(tag = "A") +
  ylim(1000, 17000) +
  annotate("text", x = 10, y = 15000, size = 5, family = "Arial", 
           label ="paste(italic(R^2), \" = 0.61 ***\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(color = "black", size = 0.5),
        plot.tag = element_text(size = 20, family = "Arial", color = "black"))

# Supplementary Fig. 5B
summary(lm(Prec_s ~ MAP_mm, data = NEON_stab.data))

Prec_s_map.plot <- ggplot(NEON_stab.data, aes(y = Prec_s, x = MAP_mm, colour = MAP_mm)) +
  geom_point(aes(size = MAT_C), alpha = 0.5) +
  scale_color_gradient(low = "orange", high = "blue") +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ x, size = 1.5, colour ="black") +
  ylab(bquote(atop(paste("Precipitation seasonality"), (paste("CV"))))) +
  xlab("MAP (mm)") +
  labs(tag = "B") +
  annotate("text", x = 1500, y = 67, size = 5, family = "Arial", 
           label ="paste(italic(R^2), \" = 0.17 *\")", parse=TRUE, hjust= 0) +
  theme_bw() +
  theme(#legend.position = "bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size = 0.5),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 16, family = "Arial", color = "black"),
        axis.text.x = element_text(size = 16, family = "Arial", color = "black"),
        axis.title.y = element_text(size = 18, family = "Arial", color = "black"),
        axis.title.x = element_text(size = 18, family = "Arial", color = "black"),
        axis.ticks.length = unit(1, "mm"),
        axis.line = element_line(color = "black", size = 0.5),
        plot.tag = element_text(size = 20, family = "Arial", color = "black"))

