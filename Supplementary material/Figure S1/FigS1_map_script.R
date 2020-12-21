############################
#### Eilat map making script
############################
## Dec 2019 by NRE, latest edit - 27/01/2020

install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", "ggspatial", "libwgeom", "sf", "rnaturalearth", "rnaturalearthdata"))

library("ggplot2")
library("ggmap")
library("sf")
library("rnaturalearth")
library("rnaturalearthdata")
library("ggspatial")
library("ggsn")
library("rgeos")

world <- ne_countries(scale = "medium", returnclass = "sf")
class(world)
world_points<- st_centroid(world)
world_points <- cbind(world, st_coordinates(st_centroid(world$geometry)))
redsea_countries<-world_points[world_points$name %in% c("Israel", "Egypt", "Saudi Arabia", "Jordan","Sudan","Eritrea", "Yemen", "Djibouti"), ]

ggplot(data = world) +
  geom_sf(color = "black", fill = "antiquewhite2") +
  xlab("Longitude") + ylab("Latitude") + 
  coord_sf(xlim = c(20, 60), ylim = c(0, 35), expand = FALSE) + 
  geom_text(data= redsea_countries,aes(x=X, y=Y, label=name), color = "black", fontface = "bold", check_overlap = FALSE,size = 4) +
  annotate(geom = "text", x = 37.8, y = 22, label = "Red Sea", angle = -60, fontface = 4, color = "grey22", size = 4) +
  annotate(geom = "text", x = 54, y = 5, label = "Indian Ocean", fontface = 4, color = "grey22", size = 3.7) +
  theme(panel.grid.major = element_blank(), panel.background = element_rect(fill = "steelblue1")) +
  geom_point(x = 34.9, y = 29.45, shape=21,fill = "red", size = 3) + theme(axis.text = element_text(colour = "black")) +  
  theme(axis.ticks.length=unit(.2,"cm")) + theme(panel.background = element_rect(colour = "black", size=1)) +
  annotation_north_arrow(which_north = "true",style = north_arrow_nautical, location = "br") +
  annotation_scale(location = "br", width_hint = 0.5, )
  
#### Map was then cleaned up using Illustrator.

