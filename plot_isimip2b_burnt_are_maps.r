graphics.off()
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/")
library("raster")
sourceAllLibs("libs/")

datFile = '../fireMIPbenchmarking/data/benchmarkData/MODIS250_q_BA_regridded0.5.nc'

dirs = paste0('../ConFIRE_ISIMIP/INFERNOonOff/Global/',  c("historic_on/"))#, "RCP6.0_on/"))
layers = list(c(1669:1728))#, 1:(17*12))

models = list.files(dirs[1])

cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')

lims = c(0, 0.1, 1, 2, 5, 10, 20, 50)

dlims = c(-10, -5, -2, -1, -0.5, 0.5, 1, 2, 5, 10)
dcols = rev(c('#b2182b','#d6604d','#f4a582','#fddbc7','#f7f7f7','#d1e5f0','#92c5de','#4393c3','#2166ac'))

diff = TRUE

openPlotBurntArea <- function(model) {
    openFile <- function(dir, layer) {
        files = list.files(paste0(dir, model), full.name = T)
        file  = files[grepl('burnt', files)]
        if (length(file) == 0) return()
        dat = brick(file)[[layer]]
        return(dat)
    }
    out = mapply(openFile, dirs, layers)
    if (length(out) == 1) out = out[[1]]
    else browser()
    out = mean(out)*12*60*60*24*30*100
    return(out)
}

#obs = mean(brick(datFile))*12*100
#sims = lapply(models, openPlotBurntArea)
#obs = raster::resample(obs, sims[[1]])
#obs[is.na(sims[[1]])] = NaN
png("figs/burnt_area_isimip2b_hist_map.png", height = 6, width = 7.2, units = 'in', res = 300)
layout(rbind(c(1, 0), c(2, 0), c(3, 4), c(5, 6), 7), heights = c(1, 0.25, 1, 1, 0.25))
par(mar = c(0, 0, 0.5, 0), oma = c(0.5, 0, 0.5, 0))

plotMap <- function(..., txt = '') {
    
    plotStandardMap(..., projection = "robinson", ylim = c(-60, 90), y_range = c(-60, 90), 
                    interior = FALSE)
    mtext(side = 3, line = -0.5, adj = 0.01, txt)
}

plotMap(obs, cols, lims, txt = 'a) Fire CCI')
legendColBar(c(0.1, 0.8), c(0.05, 0.95), cols = cols, limits = lims,
             transpose = TRUE, extend_min = FALSE, units = '%')

if (diff) {
    dsims = lapply(sims, function(sim) sim - obs)
    extend_min = TRUE
} else {
    dlims = lims
    dcols = cols
    dsims = sims
    extend_min = FALSE
}

mapply(plotMap, dsims, txt = paste0(letters[2:5], ') ', models), MoreArgs = list(dcols, dlims))
legendColBar(c(0.1, 0.8), c(0.05, 0.95), cols = dcols, limits = dlims, 
             transpose = TRUE, extend_min = extend_min, units = '%')
dev.off()

