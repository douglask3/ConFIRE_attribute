climate = list(dir = 'outputs/amazon_region/climate/',
               files = c(pr = 'precip2001-2019.nc',
                         emc = 'emc-2001-2019.nc',
                         air = 'tasMax_2001-2019.nc',
                         air = 'air2001-2019.nc',
                         soilw = 'soilw.0-10cm.gauss.2001-2019.nc'),
               cols = c("#000099", "#AA00AA", "black", "green", "cyan"))


human = list(dir = 'outputs/amazon_region/human/',
             files = c(crop = 'cropland2001-2019.nc',
                         pasture = 'pasture2001-2019.nc',
                         popden = 'population_density-PD_HYDEv3.2_2001-2018.nc'),
             cols = c("green", "brown", "black"))

vegetation = list(dir = "outputs/amazon_region/vegetation/",
             files = c(soilw_max = 'MaxOverMean_soilw.0-10cm.gauss.2001-2019.nc',
                         tree = 'treecover-2001-June2018.nc',
                         veg = 'vegcover-2001-June2018.nc'),
             cols = c("cyan", "#000099", "green"))

fireCount = list(dir = 'outputs/amazon_region/fire_counts/',
                 files = c(fireCount = 'firecount_TERRA_M__T.nc'),
             cols = c("red"))

lat = -9
lon = -65

fireMonths = 8:9

months = lapply(fireMonths, seq, 210, by = 12)
monthsByYr = lapply(1:min(sapply(months, length)),
                    function(yr) sapply(months, function(m) m[yr]))

loadDat <- function(file, dir) {
    dat = brick(paste0(dir,'/', file))
    dat = layer.apply(monthsByYr, function(i) (dat[[i]]))
    dat = dat[cellFromXY(dat, xy = cbind(lon, lat))]
    dat = dat - min(dat) + 0.00001
    dat = dat/max(dat)
    return(t(dat))
}

plotWindow <- function(info) {
    dats = lapply(info$files, loadDat, info$dir)

    plot(c(2001, 2017), c(0.0001, 1), type = 'n')
    
    mapply(function(x, col) lines(seq(2001,2017,length.out = length(x)),
           x, col = col, lwd = 2), dats, info$cols)
    lines(seq(2001,2017,length.out = length(dat_fire)), dat_fire, col = 'red', , lwd = 2, lty = 2) 
    plot(c(0,1), c(0,1), type = 'n')
    mapply(function(x, col) points(x, dat_fire, col = col, pch = 19), dats, info$cols)
}

par(mfrow = c(3,2))
dat_fire = loadDat(fireCount$files, fireCount$dir)
plotWindow(climate)
plotWindow(human)
plotWindow(vegetation)
