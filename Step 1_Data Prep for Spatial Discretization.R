### Data Prep for Time Segmentation Model ###

library(dplyr)
library(ggplot2)
library(lubridate)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(sp)
library(raster)
library(rgdal)
library(adehabitatLT)


#####################
#### Import Data ####
#####################

dat<- read.csv("gps_pos_drawdown2020_01_16JC.csv", header = T, sep = ",")

#explore
str(dat)
summary(dat)
unique(dat$id) %>% length() # number of IDs
# dat$id<- as.factor(dat$id)
names(dat)[4]<- "time"
dat$time<- as.POSIXct(strptime(dat$time, format = "%m/%d/%y %H:%M"))

#create class ltraj for SL/TA/NSD (i.e. dist, rel.angle, R2n)
dat.spdf<- dat
coordinates(dat.spdf)<- ~utmlong + utmlat
proj4string(dat.spdf)<- CRS("+init=epsg:32617")
dat.traj<- as.ltraj(xy = coordinates(dat.spdf), date = dat.spdf$time, id = dat.spdf$id)
plot(dat.traj)
dat.traj


##Assign newly calculated vars to DF
dat<- ld(dat.traj) #turn ltraj object into DF
dat<- dat[,c(11,1:10)]


#Load world map data
usa <- ne_states(country = "United States of America", returnclass = "sf")
fl<- usa %>% filter(name == "Florida") %>% st_transform(proj4string(dat.spdf))

# lakes
lakes10 <- ne_download(scale = 10, type = 'lakes', category = 'physical', returnclass = "sf")
lakes10<- sf::st_transform(lakes10, crs = "+init=epsg:32617") %>%
  sf::st_crop(xmin = min(dat$x-20000), xmax = max(dat$x+20000),
              ymin = min(dat$y-20000), ymax = max(dat$y+20000))


##########################################################
#### Create Grid to Discretize Space ####
##########################################################

# 5 km w 1 cell buffer
grid_5<- raster(extent(dat.spdf) + 10000)
res(grid_5)<- 5000
proj4string(grid_5)<- CRS("+init=epsg:32617")

grid_5[]<- 0
dat$grid.cell<- cellFromXY(grid_5, dat.spdf)

### Write to CSV for further analysis
write.csv(dat, "Snail Kite Gridded Data_TOHO.csv", row.names = F)


### Plot all points over grid

#Create grid cell borders
borders_5<- rasterToPolygons(grid_5, dissolve = F)
borders_5f<- fortify(borders_5)

#Calc points per cell
tab<- table(cellFromXY(grid_5, dat.spdf))
grid_5[as.numeric(names(tab))] <- tab
grid_5f<- as.data.frame(grid_5, xy = TRUE)
names(grid_5f)[3]<- "count"


#plot points over grid
ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_path(data = borders_5f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat, aes(x=x, y=y, color=as.factor(id)), size=0.5, alpha=0.5) +
  scale_color_viridis_d("ID", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

#plot density surface of points in grid
ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_tile(data=grid_5f, aes(x=x, y=y, fill=count)) +
  geom_path(data = borders_5f, aes(x=long, y=lat, group=group), size=0.25) +
  scale_fill_viridis_c("# of Observations", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()

#plot all tracks
ggplot() +
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_path(data = dat, aes(x=x, y=y, color = id), size = 0.5) +
  # geom_point(data = dat, aes(x=x, y=y, color = id), size = 2) +
  scale_color_viridis_d("") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()
