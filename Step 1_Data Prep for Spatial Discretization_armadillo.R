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

dat<- read.csv("three_banded_over20_for_Josh.csv", header = T, sep = ",")

#explore
str(dat)
summary(dat)
unique(dat$id) %>% length() # number of IDs
dat$date<- as.POSIXct(strptime(dat$date, format = "%d/%m/%Y %H:%M"))

#create class ltraj for SL/TA/NSD (i.e. dist, rel.angle, R2n)
dat.spdf<- dat
coordinates(dat.spdf)<- ~x + y
proj4string(dat.spdf)<- CRS("+init=epsg:32721")
dat.traj<- as.ltraj(xy = coordinates(dat.spdf), date = dat.spdf$date, id = dat.spdf$id,
                    infolocs = dat.spdf@data[,c("alt","InBurrow")])
plot(dat.traj)
dat.traj


##Assign newly calculated vars to DF
dat.traj.df<- ld(dat.traj) #turn ltraj object into DF
dat.new<- dat.traj.df[,c(11,1:10,13:14)]



#north region
dat.N<- dat[dat$id %in% c('tm14','tm15','tm21','tm30'),]

ggplot() +
  geom_path(data = dat.N, aes(x = x, y = y, color = id)) +
  scale_color_viridis_d("ID") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  coord_equal()

#south region
dat.S<- dat[dat$id %in% c('tm2','tm3','tm4','tm5','tm8'),]

ggplot() +
  geom_path(data = dat.S, aes(x = x, y = y, color = id)) +
  scale_color_viridis_d("ID") +
  labs(x = "Easting", y = "Northing") +
  theme_bw() +
  coord_equal()

##########################################################
#### Create Grid to Discretize Space ####
##########################################################

### Northern group

dat.N.spdf<- dat.spdf[dat.spdf@coords[,2] > 8210000,]

# 200 m w 1 cell buffer
grid_200N<- raster(extent(dat.N.spdf) + 400)
res(grid_200N)<- 200
proj4string(grid_200N)<- CRS("+init=epsg:32721")

grid_200N[]<- 0
dat.N$grid.cell<- cellFromXY(grid_200N, dat.N.spdf)


## Plot all points over grid

#Create grid cell borders
borders_200N<- rasterToPolygons(grid_200N, dissolve = F)
borders_200Nf<- fortify(borders_200N)

#Calc points per cell
tab<- table(cellFromXY(grid_200N, dat.N.spdf))
grid_200N[as.numeric(names(tab))] <- tab
grid_200Nf<- as.data.frame(grid_200N, xy = TRUE)
names(grid_200Nf)[3]<- "count"


#plot points over grid
ggplot() +
  coord_sf(xlim = c(min(dat.N$x-250), max(dat.N$x+250)),
           ylim = c(min(dat.N$y-250), max(dat.N$y+250)), expand = FALSE) +
  geom_path(data = borders_200Nf, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat.N, aes(x=x, y=y, color=as.factor(id)), size=0.5, alpha=0.5) +
  scale_color_viridis_d("ID", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())

#plot density surface of points in grid
ggplot() +
  coord_sf(xlim = c(min(dat.N$x-250), max(dat.N$x+250)),
           ylim = c(min(dat.N$y-250), max(dat.N$y+250)), expand = FALSE) +
  geom_tile(data=grid_200Nf, aes(x=x, y=y, fill=count)) +
  geom_path(data = borders_200Nf, aes(x=long, y=lat, group=group), size=0.25) +
  scale_fill_viridis_c("# of Observations", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())






### Southern group

dat.S.spdf<- dat.spdf[dat.spdf@coords[,2] < 8210000,]

# 100 m w 1 cell buffer
grid_100S<- raster(extent(dat.S.spdf) + 200)
res(grid_100S)<- 100
proj4string(grid_100S)<- CRS("+init=epsg:32721")

grid_100S[]<- 0
dat.S$grid.cell<- cellFromXY(grid_100S, dat.S.spdf)


## Plot all points over grid

#Create grid cell borders
borders_100S<- rasterToPolygons(grid_100S, dissolve = F)
borders_100Sf<- fortify(borders_100S)

#Calc points per cell
tab<- table(cellFromXY(grid_100S, dat.S.spdf))
grid_100S[as.numeric(names(tab))] <- tab
grid_100Sf<- as.data.frame(grid_100S, xy = TRUE)
names(grid_100Sf)[3]<- "count"


#plot points over grid
ggplot() +
  coord_sf(xlim = c(min(dat.S$x-250), max(dat.S$x+250)),
           ylim = c(min(dat.S$y-250), max(dat.S$y+250)), expand = FALSE) +
  geom_path(data = borders_100Sf, aes(x=long, y=lat, group=group), size=0.25) +
  geom_point(data = dat.S, aes(x=x, y=y, color=as.factor(id)), size=0.5, alpha=0.5) +
  scale_color_viridis_d("ID", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())

#plot density surface of points in grid
ggplot() +
  coord_sf(xlim = c(min(dat.S$x-250), max(dat.S$x+250)),
           ylim = c(min(dat.S$y-250), max(dat.S$y+250)), expand = FALSE) +
  geom_tile(data=grid_100Sf, aes(x=x, y=y, fill=count)) +
  geom_path(data = borders_100Sf, aes(x=long, y=lat, group=group), size=0.25) +
  scale_fill_viridis_c("# of Observations", alpha = 0.6) +
  labs(x = "Longitude", y = "Latitude") +
  theme_bw() +
  theme(panel.grid = element_blank())







### Write to CSV for further analysis
write.csv(dat.N, "Armadillo Gridded Data_N.csv", row.names = F)
write.csv(dat.S, "Armadillo Gridded Data_S.csv", row.names = F)
