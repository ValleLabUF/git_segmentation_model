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

#rivers
rivers10 <- ne_download(scale = 10, type = 'rivers_lake_centerlines', category = 'physical',
                        returnclass = "sf")
rivers10<- sf::st_transform(rivers10, crs = "+init=epsg:32643") %>%
  sf::st_crop(xmin = min(dat$x-11000), xmax = max(dat$x+11000),
              ymin = min(dat$y-11000), ymax = max(dat$y+11000))


##########################################################
#### Create Grid to Discretize Space ####
##########################################################

res.function.optim=function(param, dat.spdf, crs){
  quant = param[1]
  
  #Set resolution as some quantile of SL
  res<- quantile(sqrt(dat.spdf$R2n), quant, na.rm = T)
  
  
  # create grid
  grid<- raster(extent(dat.spdf) + (2*res))
  res(grid)<- res
  proj4string(grid)<- crs
  grid[]<- 0
  
  time.series<- list()
  for(i in 1:length(unique(dat.spdf$id))) {
    
    time.series[[i]]<- cellFromXY(grid, subset(dat.spdf, id == unique(dat.spdf$id)[i]))
  }
  
  
  ts.length<- lapply(time.series, rle) %>% 
    lapply(., '[[', 1) %>% 
    unlist() %>% 
    mean()
  
  n.length<- lapply(time.series, rle) %>% 
    lapply(., '[[', 1) %>%
    lapply(., function(x) length(which(x > 10))) %>% 
    unlist() %>% 
    mean()
  
  -ts.length * n.length
}

snwlpd.test<- optim(par = 0.5, fn = res.function.optim, dat.spdf = dat.spdf, lower = 0.1,
                    upper = 0.999, crs = CRS("+init=epsg:32643"), method = "Brent")



quant<- snwlpd.test$par


#Set resolution as 99.5% quantile of SL
res<- quantile(sqrt(dat$R2n), quant, na.rm = T)  #4.5 km


# 4.5 km w 1 cell buffer
grid<- raster(extent(dat.spdf) + (2*res))
res(grid)<- res
proj4string(grid)<- CRS("+init=epsg:32643")

grid[]<- 0
dat$grid.cell<- cellFromXY(grid, dat.spdf)

### Write to CSV for further analysis
write.csv(dat, "Snow Leopard Gridded Data.csv", row.names = F)


### Plot all points over grid

#Create grid cell borders
borders<- rasterToPolygons(grid, dissolve = F)
borders_f<- fortify(borders)

#Calc points per cell
tab<- table(cellFromXY(grid, dat.spdf))
grid[as.numeric(names(tab))] <- tab
grid_f<- as.data.frame(grid, xy = TRUE)
names(grid_f)[3]<- "count"


#plot points over grid
ggplot() +
  geom_sf(data = afg) +
  geom_sf(data = rivers10, color = "lightblue", alpha = 0.65, lwd = 5) +
  coord_sf(xlim = c(min(dat$x-10000), max(dat$x+10000)),
           ylim = c(min(dat$y-10000), max(dat$y+10000)), expand = FALSE) +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat, aes(x=x, y=y, fill=as.factor(id)), pch = 21, size=1.5, alpha=0.5) +
  scale_fill_viridis_d("ID", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())

#plot density surface of points in grid
ggplot() +
  geom_sf(data = afg) +
  geom_sf(data = rivers10, color = "lightblue", alpha = 0.65, lwd = 5) +
  coord_sf(xlim = c(min(dat$x-10000), max(dat$x+10000)),
           ylim = c(min(dat$y-10000), max(dat$y+10000)), expand = FALSE) +
  geom_tile(data=grid_f, aes(x=x, y=y, fill=count)) +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
  scale_fill_viridis_c("# of Observations", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())

#plot all tracks
ggplot() +
  geom_sf(data = afg) +
  geom_sf(data = rivers10, color = "lightblue", alpha = 0.65, lwd = 5) +
  coord_sf(xlim = c(min(dat$x-10000), max(dat$x+10000)),
           ylim = c(min(dat$y-10000), max(dat$y+10000)), expand = FALSE) +
  geom_path(data = dat, aes(x=x, y=y, color = id), size = 0.5) +
  scale_color_viridis_d("") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())
