graphics.off()
library(raster)
source("libs/plotStandardMap.r")
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterPlotFunctions/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")

cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a',
            '#e31a1c','#bd0026','#800026')

likiDir = "outputs/sampled_posterior_ConFire_ISIMIP_solutions/attempt3-full/"
likiFiles = paste0(c("historic", "RCP2.6", "RCP6.0"), "/fullPost.nc")

ObsFile = "/prj/ukesm/doukel/LimFIRE/outputs/longTermRecord_modis--Max_monthly_burnt_areaJan_2001-Dec_2019.nc"

seamaskFile = "data/seamask.nc"

pnts = list("Artic fires 2019" = c(120, 68), "Amazon fires 2019" = c(-53, -10.5),
            "SE Aus 2019/2020" = c(147.5, -35))
pnts = list("Artic fires 2019" = c(120, 68), "Amazon fires 2019" = c(-53, -10.5),
            "SE Aus 2019/2020" = c(151, -30))#, "Coastal Argentina 2011" = c( -60,  -37.5))

obs = raster(ObsFile)/4

openLiki <- function(file) 
    lapply(paste0(likiDirs, '/', file), brick)

likiDirs = list.dirs(likiDir, recursive = FALSE)[c(1, 2)]
likis = lapply(likiFiles, openLiki)

seamask = raster(seamaskFile)

obs = raster::resample(obs, likis[[1]][[1]])
obs[seamask == 0] = NaN
plotStandardMap(obs, limits = c(0, 1, 2, 5, 10, 20, 50, 100), col = cols, y_range = c(-60, 90))
lapply(pnts, function(i) points(i[1], i[2]))

px = names(likis[[1]][[1]])

convertLayer2Num <- function(x) {
    x0 = x
    x = substr(x, 2, nchar(x))
    
    if (substr(x, 1, 1) == '.') x = - as.numeric(substr(x, 2, nchar(x)))
    else x = as.numeric(x)
    #logistic(x)
}
px = sapply(px, convertLayer2Num)
forPnt <- function(pnt, nm) {
    xy = c(colFromX(obs, pnt[1]), rowFromY(obs, pnt[2]))
    
    addPoly <- function(liki, col = NULL) {
        
        py = lapply(liki, function(i) i[xy[2], xy[1]])
        py = py0 = lapply(py, function(i) i / sum(i))
        py = do.call('*', py)^(1/length(py))
        #lapply(py, function(y) polygon(px, y, col = make.transparent(col, 0.9), border = NA))
  
        if (is.null(col))     
            return(unlist(py))
        else 
            polygon(px, py, col = make.transparent(col, 0.95), border = NA)

        return(py)
    }
    py = mapply(addPoly, likis)
    test = which(apply(py, 1, sum)>0)
    test = c(max(1, test[1]-1), test, min(length(py), tail(test, 1)+1))
    xr = range(px[test])
    #xr = range(px)
    xr = logit(c(0.1, 30)/100)
    plot(xr, c(0, max(unlist(py))),
         type = 'n', axes = FALSE, xlab = '', ylab = '')
    for (i in 1:8) mapply(addPoly, rev(likis), c("red", "blue", "black"))
    #polygon(px, py[,2], col = make.transparent('blue', 0.9), border = 'blue')
    #polygon(px, py[,3], col = make.transparent('red', 0.9), border = 'red')
    #polygon(px, py[,1], col = make.transparent(, 0.9))
    
   
    labels = c(0, signif(logistic(seq(xr[1], xr[2], length.out = 5)), 1), 1)
    at = logit(labels)
    labels = labels *100
    at[1] = min(px); at[length(at)] = max(px)
    axis(1, at = at, labels = labels)
    
    ii = xy[1] + -1:1
    jj = xy[2] + -1:1
    #browser()
    pO = p0 = obs[jj, ii]
    pO = logit(mean(pO)/100) - 0.0
    lines(c(pO, pO), c(0, max(py)*0.67), lty = 2)
    pc =  apply(py, 2, function(i) sum(i[px > pO])/sum(i))
    if (pc[1] < 0.01) pc[1] = 0.01
    pc[2:3] = pc[2:3]/pc[1]
    pc[1] = pc[1] * 100
    pc = round(pc, 2)
    if (pc[1] <= 1) {
        pc[2:3] =  round((100/pc[2:3])-1)
        pc =  as.character(pc)       
        pc[1] = '< 1'   
        
        mtext(side = 3, line = -3, adj = 1, col = 'blue', 
              paste0('occurs one in ', pc[2], '  years by 2100 under RCP2.6'))
        mtext(side = 3, line = -4.5, adj = 1, col = 'red',  
              paste0('occurs one in ', pc[3], '  years by 2100 under RCP6.0s'))
    
    } else {
        mtext(side = 3, line = -3, adj = 1, col = 'blue', 
              paste0(pc[2], ' times more likely by 2100 under RCP2.6'))
        mtext(side = 3, line = -4.5, adj = 1, col = 'red',  
              paste0(pc[3], ' times more likely by 2100 under RCP6.0'))
    }
        
    #browser()
    #pc =
    #pc[is.infinite(pc)] = '100+'
    mtext(side = 3, line = -1.5, adj = 1,  paste0(nm, ' at ', pc[1], '% likihood'))

    browser()
    
}

png("figs/likihoodEventCurves.png", height = 6, width = 5, units = 'in', res = 300)
    par(mfrow = c(3, 2), mar = rep(1, 4), oma = c(2, 2, 0, 0))
    mapply(forPnt,pnts, names(pnts))
    mtext('Burnt area (%)', side = 1, line = 2, font = 2)
    mtext('Prob.', side = 2, line = 0, outer = TRUE, font = 2)
dev.off()
