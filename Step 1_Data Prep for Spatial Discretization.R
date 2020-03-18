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


#########################################
#### Create Grid to Discretize Space ####
#########################################


res.function=function(quant, dat.spdf, crs){
  
  grid.res<- matrix(NA, length(quant), 3)
  colnames(grid.res)<- c("res","n.cells", "n.length")
  grid.res[,1]<- quantile(sqrt(dat.spdf$R2n), quant, na.rm = T)
  
  for (i in 1:length(quant)) {
    #Set resolution as some quantile of displacement
    res<- quantile(sqrt(dat.spdf$R2n), quant[i], na.rm = T)
    
    # create grid
    grid<- raster(extent(dat.spdf) + (2*res))
    res(grid)<- res
    proj4string(grid)<- crs
    grid[]<- 0
    
    
    time.series<- list()
    for(j in 1:length(unique(dat.spdf$id))) {
      
      time.series[[j]]<- cellFromXY(grid, subset(dat.spdf, id == unique(dat.spdf$id)[j]))
    }
    
    ind.occup.cells<- sum(sapply(time.series, function(x) length(unique(x))) > 2)
    all.occup.cells<- unlist(time.series) %>% unique() %>% length()
    
    n.cells<- time.series %>% 
      lapply(., function(x) length(unique(x))) %>% 
      unlist() %>% 
      median()
    
    n.length<- lapply(time.series, rle) %>%
      lapply(., '[[', 1) %>%
      lapply(., median) %>% 
      unlist() %>% 
      median()
    
     
    grid.res[i,2]<- n.cells
    grid.res[i,3]<- n.length
  }
  
  grid.res
}

dat.spdf@data<- dat

quant<- seq(0.2, 0.5, length.out = 120)
quant.range<- res.function(quant = quant, dat.spdf = dat.spdf, crs = proj4string(dat.spdf))
quant.range[,2]<- quant.range[,2]/max(quant.range[,2])
quant.range[,3]<- quant.range[,3]/max(quant.range[,3])

quant.range<- quant.range %>%
  data.frame() %>%
  mutate(dist = sqrt((n.cells-0)^2 + (n.length-0)^2))

ggplot() +
  geom_point(data = quant.range, aes(x=n.cells, y=n.length, color = dist), size = 2,
             alpha = 0.8) +
  coord_equal() +
  scale_color_viridis_c(direction = -1) +
  theme_bw()

quant.range2<- quant.range %>% top_n(n=-10, wt=dist)  #12.7 - 28.3 km



res<- 15000  #min value from quant fun

# 15 km w 1 cell buffer
grid<- raster(extent(dat.spdf) + 2*res)
res(grid)<- res
proj4string(grid)<- CRS("+init=epsg:32617")

grid[]<- 0
dat$grid.cell<- cellFromXY(grid, dat.spdf)

### Write to CSV for further analysis
write.csv(dat, "Snail Kite Gridded Data_TOHO.csv", row.names = F)


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
  geom_sf(data = fl) +
  geom_sf(data = lakes10, fill = "lightblue", alpha = 0.65) +
  coord_sf(xlim = c(min(dat$x-20000), max(dat$x+20000)),
           ylim = c(min(dat$y-20000), max(dat$y+20000)), expand = FALSE) +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
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
  geom_tile(data=grid_f, aes(x=x, y=y, fill=count)) +
  geom_path(data = borders_f, aes(x=long, y=lat, group=group), size=0.25) +
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
  scale_color_viridis_d("") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw()
