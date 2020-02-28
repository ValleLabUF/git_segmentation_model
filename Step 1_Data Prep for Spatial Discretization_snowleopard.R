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

dat<- read.csv("Modified Snow Leopard Data.csv", header = T, sep = ",")

#explore
str(dat)
summary(dat)
unique(dat$id) %>% length() # number of IDs
dat$date<- dat$date %>% as_datetime()

#create SPDF object
dat.spdf<- dat
coordinates(dat.spdf)<- ~x + y
proj4string(dat.spdf)<- CRS("+init=epsg:32643")


#Load world map data
afg <- ne_states(country = "Afghanistan", returnclass = "sf") %>%
  st_transform(proj4string(dat.spdf))
# fl<- usa %>% filter(name == "Florida") %>% st_transform(proj4string(dat.spdf))

#rivers
rivers10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical',
                        returnclass = "sf")
rivers10<- sf::st_transform(rivers10, crs = "+init=epsg:32643") %>%
  sf::st_crop(xmin = min(dat$x-5500), xmax = max(dat$x+5500),
              ymin = min(dat$y-5500), ymax = max(dat$y+5500))


##########################################################
#### Create Grid to Discretize Space ####
##########################################################

res.function.optim=function(param, dat.spdf, crs){
  quant = param[1]
  
  #Set resolution as some quantile of SL
  res<- quantile(dat.spdf$dist, quant, na.rm = T)  #3.8 km
  
  
  # create grid
  grid<- raster(extent(dat.spdf) + (2*res))
  res(grid)<- res
  proj4string(grid)<- crs
  grid[]<- 0
  
  time.series<- list()
  for(i in 1:length(unique(dat.spdf$id))) {
    
    time.series[[i]]<- cellFromXY(grid, subset(dat.spdf, id == unique(dat.spdf$id)[i]))
  }
  
  
  lapply(time.series, function(x) mean(diff(x)^2)) %>%
    unlist() %>%
    mean()  #Mean difference between grid cell numbers averaged across all IDs
}

snwlpd.test<- optim(par = 0.95, fn = res.function.optim, dat.spdf = dat.spdf, lower = 0.9,
                    upper = 0.99, crs = CRS("+init=epsg:32643"), method = "Brent")



quant<- snwlpd.test$par


#Set resolution as 95% quantile of SL
res<- quantile(dat$dist, quant, na.rm = T)  #3.5 km


# 3.5 km w 1 cell buffer
grid_5<- raster(extent(dat.spdf) + (2*res))
res(grid_5)<- res
proj4string(grid_5)<- CRS("+init=epsg:32643")

grid_5[]<- 0
dat$grid.cell<- cellFromXY(grid_5, dat.spdf)

### Write to CSV for further analysis
write.csv(dat, "Snow Leopard Gridded Data.csv", row.names = F)


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
  geom_sf(data = afg) +
  geom_sf(data = rivers10, color = "lightblue", alpha = 0.65, lwd = 5) +
  coord_sf(xlim = c(min(dat$x-5000), max(dat$x+5000)),
           ylim = c(min(dat$y-5000), max(dat$y+5000)), expand = FALSE) +
  geom_path(data = borders_5f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat, aes(x=x, y=y, fill=as.factor(id)), pch = 21, size=1.5, alpha=0.5) +
  scale_fill_viridis_d("ID", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())

#plot density surface of points in grid
ggplot() +
  geom_sf(data = afg) +
  geom_sf(data = rivers10, color = "lightblue", alpha = 0.65, lwd = 5) +
  coord_sf(xlim = c(min(dat$x-5000), max(dat$x+5000)),
           ylim = c(min(dat$y-5000), max(dat$y+5000)), expand = FALSE) +
  geom_tile(data=grid_5f, aes(x=x, y=y, fill=count)) +
  geom_path(data = borders_5f, aes(x=long, y=lat, group=group), size=0.25) +
  scale_fill_viridis_c("# of Observations", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())

#plot all tracks
ggplot() +
  geom_sf(data = afg) +
  geom_sf(data = rivers10, color = "lightblue", alpha = 0.65, lwd = 5) +
  coord_sf(xlim = c(min(dat$x-5000), max(dat$x+5000)),
           ylim = c(min(dat$y-5000), max(dat$y+5000)), expand = FALSE) +
  geom_path(data = dat, aes(x=x, y=y, color = id), size = 0.5) +
  # geom_point(data = dat, aes(x=x, y=y, color = id), size = 2) +
  scale_color_viridis_d("") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())
