##### Basic data prep

# clean working environment
rm(list = ls(all.names = T))
gc()

## prevent encoding error
Sys.getlocale()
Sys.setlocale("LC_CTYPE", ".1251")
Sys.getlocale()

# load packages
library(terra)
library(dplyr)


##### prep environment data  ----------

# buffer
buff <- terra::vect('buffer/buffer_300km.shp')

# climate 
clim <- rast(list.files(path = 'E:/env layers/worldclim', pattern = '.tif$', full.names = T))
plot(clim[[1]])

# elev
elev <- rast('E:/env layers/elev_worldclim/wc2.1_30s_elev.tif')

# stack, crop, mask
envs <- c(clim, elev)
envs <- crop(envs, ext(buff))
envs <- mask(envs, buff)
plot(envs[[1]])

# rename files
names(envs) = c('bio1','bio10','bio11', 'bio12','bio13','bio14','bio15','bio16','bio17',
                'bio18','bio19','bio2','bio3','bio4','bio5','bio6','bio7','bio8','bio9','elev')

# export files
for (i in 1:nlyr(envs)) {
  r <- envs[[i]]
  name <- paste0('envs/', names(envs)[[i]], '.tif')
  writeRaster(r, filename = name, overwrite = T)
}


##### sample background points ----------
# random sampling
set.seed(333)

bg <- dismo::randomPoints(mask = raster::raster(envs[[1]]), n = 10000) %>% as.data.frame()
colnames(bg) = c('long', 'lat')
head(bg)

# export
write.csv(bg, 'bg/bg_10000.csv')


##### run correlation test ----------
ntbox::run_ntbox()

# Spearman's test |r| > 0.7 cutoff ==  bio1 bio12 bio14 bio2 bio3 


