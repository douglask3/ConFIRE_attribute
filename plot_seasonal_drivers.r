graphics.off()
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/")
library("raster")

gswp3_path = "/hpc/data/d00/hadea/isimip3a/InputData/climate/atmosphere/obsclim/GSWP3-W5E5/gswp3-w5e5_obsclimfill_pr_global_daily_"

gswp3_years = c("1971_1980", "1981_1990", "1991_2000", "2001_2010", "2011_2019")
obs_years = 1979:2013

openGSWP3 <- function(yr) {
    print(yr)
    dat = brick(paste0(gswp3_path, yr, '.nc'))
    date = names(dat)
    test = sapply(date, function(i) any(substr(i, 2, 5) == obs_years))    
    return(dat[[which(test)]])
}

#gswp3 = layer.apply(gswp3_years, openGSWP3)
nlperfile = 100

sumLayers <- function(i) {
    print(i)
    tfile = paste0("temp/gswp3_", i, ".nc")
    if (file.exists(tfile)) return(raster(tfile))
    index = i:min(nlayers(gswp3), i+nlperfile-1)
    out = sum(gswp3[[index]])
    writeRaster(out, tfile, overwrite = TRUE)

}

#gswp3i = layer.apply(seq(1, nlayers(gswp3), by = nlperfile), sumLayers)
#gswp3j = sum(gswp3i)/nlayers(gswp3)

#gswp3j = gswp3j*60*60*24


dir = '../ConFIRE_ISIMIP/INFERNOonOff/Global/'

expDir = list("No fire" = c("historic_off/", "RCP6.0_off/"), 
              "With fire" = c("historic_on/", "RCP6.0_on/"))

layers = list(1429:1728, 1:48)
Models = c("GFDL-ESM2M", "HADGEM2-ES", "IPSL-CM5A-LR", "MIROC5")

vars = c("Pr" = "precip", "SW" = "sw_down", "T" = "temp",
         "Soil Moisture" = "smc", "NPP" = "npp")
scales = c(60*60*24*30, 1, 1, 1, 1)
shifts = c(0, 0, -273.15, 0, 0)

units = c("mm", "W ~m-2~", "~DEG~C", "???", "gC ~m-2~")
sides = c(1, 3, 3, 1, 1)

mask = raster("data/basins.nc")==1
mask[is.na(mask)] = 0

climatology <- function(r) {
    climMean <- function(mn) mean(r[[seq(mn, nlayers(r), by = 12)]])
    out = layer.apply(1:12, climMean)
    out
}

rAv <- function(r, mask) {
    areaR =  raster::area(r, na.rm = TRUE)
    sum(r[mask] * areaR[mask], na.rm = TRUE)/sum(areaR[mask], na.rm = TRUE)
}

openExp <- function(exp, name, var, scale, shift, outTS = TRUE, outClim = TRUE) {
    
    openMods <- function(mod) {
        files = paste0(dir, exp, mod, '/', var, '.nc')
        tfile = paste("temp/climatology", name, mod, var, ".nc", sep = '-')
        print(tfile)
        if (file.exists(tfile)) dat = brick(tfile) else {
            dat = mapply(function(f, l) brick(f)[[l]], files, layers)
            dat = addLayer(dat[[1]], dat[[2]])
            dat = climatology(dat)
            dat = writeRaster(dat, file = tfile, overwrite = TRUE)
        }
            
        if (outTS) dat = unlist(layer.apply(dat, rAv, mask) ) else dat = mean(dat)
        dat
    }
    tfile = paste('temp/climateology-ts', name, var, outTS, '.Rd', sep = '-')
    if (file.exists(tfile)) load(tfile) else {
        dats = lapply(Models, openMods)
        save(dats, file = tfile)
    }
        
    dats = lapply(dats, function(i) scale * (shift + i))
}

prDat = mapply(openExp, expDir, names(expDir), 
                 MoreArgs = list(vars[1], 60*60*24, 0, FALSE))[,1]


png("figs/precip.png", height = 7*1.3, width = 6*1.3, units = 'in', res = 300)
layout(cbind(1:6, c(0, 7:11)), height = c(1, 1, 1, 1,1,0.3))
par(mar = rep(0, 4), oma = c(2, 1, 0, 0))
cols = c('#ffffd9','#edf8b1','#c7e9b4','#7fcdbb','#41b6c4','#1d91c0','#225ea8','#253494','#081d58')
dcols = c('#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e')
dlimits = c(-1, -0.5, -0.2, -0.1, -0.05, -0.02, -0.01, 0.01, 0.02, 0.05, 0.1, 0.20, 0.5, 1)
limits = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
gswp3j = raster::resample(gswp3j, prDat[[1]])
gswp3j[is.na(prDat[[1]])] =  NaN

plotMap <- function(r, txt, cols, limits) {
    plotStandardMap(r, cols, limits, ylim = c(-60, 90))
    mtext(txt, side = 2, adj = 0.1, line = -1)
}
plotMap(gswp3j, "GSWP3", cols, limits)

mapply(plotMap, prDat, Models, MoreArgs = list(cols, limits))
StandardLegend(cols, limits, gswp3j)

mapply(plotMap, lapply(prDat, '-', gswp3j), rep('', length(Models)), 
        MoreArgs = list(dcols, dlimits))
StandardLegend(dcols, dlimits, gswp3j, extend_min = TRUE, oneSideLabels = FALSE)

dev.off()

plotVar <- function(var, name, unit, side, scale, shift, cols = c('blue', 'red')) {
   #vars, names(vars), units, sides, scales, shifts 
    dat = mapply(openExp, expDir, names(expDir), 
                 MoreArgs = list(var, scale, shift))
    
    yrange = range(unlist(dat))
    plot(c(0, 12), yrange, xlab = '', ylab = '', xaxt = 'n', type = 'n')
    polygon(c(-9E9, -9E9, 9E9, 9E9), c(-9E9, 9E9, 9E9, -9E9), col = "white")
    axis(1, at = 0:12, col = 'white', 
         labels = c('D', 'J','F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
    axis(1, at = -0.5:12.5, labels = rep('', 14))
    plotExp <- function(dat, col) {
        col =  make.transparent(col, 0.7)
        plotMod <- function(y, lty) 
            lines(-1:14, y[c(11:12, 1:12, 1:2)], col = col, lty = lty, cex = 1.5)
        
        mapply(plotMod, dat, 1:4)
    }
    for (nn in 1:10) {
        plotExp(dat[,1], cols[1])
        plotExp(dat[,2], cols[2])
        plotExp(dat[,2], cols[2])
        plotExp(dat[,1], cols[1])
    }
    mtext(side = side, adj = 0.1, line = -1.5, name)
    mtext.units(side = 2, line = 2, unit)
    return(dat)
}

nplots = length(vars)+1
nrows = ceiling(sqrt(nplots))
ncols = ceiling(nplots/nrows)
png("figs/Amazon_seasonal.png", height = 7, width = 7, units = 'in', res = 300)
par(mfrow = c(nrows, ncols), mar = c(0.25, 3.5, 0.25, 0.5), oma = c(2.5, 0.7, 1, 0))
outs = mapply(plotVar, vars, names(vars), units, sides, scales, shifts)

plot.new()
par(mar = c(0, 5, 0, 5))
legend('right', names(expDir), lty = 1, col = c('blue', 'red'), bty = 'n', lwd = 1.5)
legend('left', Models, lty = 1:4, bty = 'n', lwd = 1.5)
dev.off()

png("figs/Amazon_seasonal_pr.png", height = 5, width = 7, units = 'in', res = 300)
#par(mfrow = c(2, 1), mar = c(0.25, 3.5, 0.25, 0.5), oma = c(2.5, 0.7, 1, 0))
outs = plotVar( vars[1], '', units[1], sides[1], scales[1], shifts[1], 
               cols = c('black', 'black'))

#plot.new()
#par(mar = c(0, 5, , 5))
#legend('right', names(expDir), lty = 1, col = c('blue', 'red'), bty = 'n', lwd = 1.5)
legend('bottomleft', Models, lty = 1:4, bty = 'n', lwd = 1.5)
dev.off()
