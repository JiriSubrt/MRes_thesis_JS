# Map for svalbard and sites
# Jiri Subrt
# 27/05/24

# Install and load libraries
#install.packages("ggspatial")
library(ggOceanMaps)
library(ggspatial)
library(patchwork)

# Load data
map_data <- read.csv("data/map_data.csv")

# Map 1: Arctic circumpolar
arctic.map <- basemap(limits = 65, land.col = "#c9c7c1")

# Map 2: Svalbard full
svalbard.map <- basemap(limits = c(11, 35, 76, 81), rotate = TRUE, land.col = "#c9c7c1") +
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") +
  #ggspatial::geom_spatial_point(data = map_data, aes(x = long, y = lat), color = "red") +
  #annotate(geom = "rect", xmax = 13, xmin = 10, ymax = 79.20, ymin = 78.7, colour = "red", fill = NA)
  theme(axis.title = element_blank())

# Map 3: Kongsfjord zoom-in
kongsfjord_map <- basemap(limits = c(10, 13, 78.7, 79.20), rotate = TRUE, land.col = "#c9c7c1") +
  ggspatial::annotation_scale(location = "br") + 
  ggspatial::annotation_north_arrow(location = "tr", which_north = "true") +
  ggspatial::geom_spatial_point(data = map_data, aes(x = long, y = lat), color = "red") +
  annotate("text", x = map_data$long, y = map_data$lat, label = map_data$location, size = 3) +
  #annotate(geom = "rect", xmax = 13, xmin = 10, ymax = 79.20, ymin = 78.7, colour = "red", fill = NA)
  theme(axis.title = element_blank())
  
# All three maps together
# I added red rectangle in the middle map manually in PowerPoint
map_together <- arctic.map +svalbard.map + kongsfjord_map +
  plot_layout(nrow = 1, heights = c(10, 10))
  
# Save map
ggsave("outcomes/svalbard_map.png", plot = map_together, width = 10, height = 8, dpi = 300)

