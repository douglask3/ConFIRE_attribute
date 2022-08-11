library(raster)

eg = raster('../ConFIRE_ISIMIP/INFERNOonOff/Global/historic_off/GFDL-ESM2M/trees.nc')
extent = extent(eg)

eg = raster::extend(eg, extent(c(-180, 180, -90, 90)))
dat = as.matrix(read.table("/home/h06/hadng/TRIP0p5/rivnum05.asc"))
eg[] = dat
eg = raster::crop(eg, extent)
eg[eg == -999] = NaN

writeRaster(eg, file = 'data/basins.nc')
