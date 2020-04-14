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
levels(dat$id)[4:5]<- "Pari"

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
  sf::st_crop(xmin = min(dat$x-21000), xmax = max(dat$x+21000),
              ymin = min(dat$y-21000), ymax = max(dat$y+21000))


##########################################################
#### Create Grid to Discretize Space ####
##########################################################


res.function=function(quant, dat.spdf, crs){
  
  # grid.list<- list()
  
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
    
    # grid.list[[i]]<- ifelse(all.occup.cells < 2000 &
    #                           (all.occup.cells/ncell(grid)) >= 0.1 &
    #                           (ind.occup.cells/length(time.series)) >= 0.75 &
    #                           (ind.n.length/length(time.series)) >= 0.75,
    #                         quant[i],
    #                         NA)
    
    grid.res[i,2]<- n.cells
    grid.res[i,3]<- n.length
  }
  
  # ind<- unlist(grid.list)
  # c(min(ind, na.rm = T), max(ind, na.rm = T))
  grid.res
}


quant<- seq(0.01, 1.0, length.out = 100)
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

quant.range2<- quant.range %>% top_n(n=-10, wt=dist)  #9.2 - 21.2 km







#Set resolution
res<- 7000  #min value of range
buffer<- 2*res

# 7 km w 1 cell buffer
grid<- raster(extent(dat.spdf) + buffer)
res(grid)<- res
proj4string(grid)<- CRS("+init=epsg:32643")

grid[]<- 0
dat$grid.cell<- cellFromXY(grid, dat[,c("x","y")])

### Write to CSV for further analysis
# write.csv(dat, "Snow Leopard Gridded Data.csv", row.names = F)


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
  coord_sf(xlim = c(min(dat$x-buffer), max(dat$x+buffer)),
           ylim = c(min(dat$y-buffer), max(dat$y+buffer)), expand = FALSE) +
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
  coord_sf(xlim = c(min(dat$x-buffer), max(dat$x+buffer)),
           ylim = c(min(dat$y-buffer), max(dat$y+buffer)), expand = FALSE) +
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
  coord_sf(xlim = c(min(dat$x-buffer), max(dat$x+buffer)),
           ylim = c(min(dat$y-buffer), max(dat$y+buffer)), expand = FALSE) +
  geom_path(data = dat, aes(x=x, y=y, color = id), size = 0.5) +
  scale_color_viridis_d("") +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())
