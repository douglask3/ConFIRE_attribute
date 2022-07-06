source("libs/filename.noPath.r")
library(ncdf4)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")

sourceAllLibs("libs/")
options(error=recover)
graphics.off()

dir = '../ConFIRE_ISIMIP/INFERNOonOff/Global/'

historicID = "historic"
futuresID = c("RCP2.6", "RCP6.0")

experiments = c("_on", "_off")
models = c("GFDL-ESM2M", "HADGEM2-ES", "MIROC5", "IPSL-CM5A-LR")

variables = c("trees", 'cveg', 'csoil') #, '
FUNs = list(function(x) logit(0.00005 + 0.999 * x), 
            function(x) log(x + 0.00005), function(x) log(x + 0.00005))

nsigYrs = 20


sigTS <- function(model, variable, FUN) {
    
    #print(tempFile)
    #if (file.exists(tempFile)) {
    #    load(tempFile)
    #    return(out)
    #}#
    openDat <- function(period, experiment) 
        brick(paste0(dir, '/', period, experiment, '/', model, '/', variable, '.nc'))

    openRatio <- function(sl, period) {
        tempFile = paste('temp/fireImpact-cells', model, variable, period, 
                         sl, ".nc", sep = '-')
        if (file.exists(tempFile)) return(raster(tempFile))
        print(tempFile)
        index = sl:(sl+12-1)
        on = sum(openDat(period, experiments[1])[[index]])
        off = sum(openDat(period, experiments[2])[[index]])

        out = on/max(addLayer(on, off))

        out = writeRaster(out, file = tempFile, overwrite = TRUE)
        return(out)
    }
    nHist = nlayers(openDat(historicID, experiments[1]))
    nFutr = nlayers(openDat(futuresID[1], experiments[1]))
    
    
    openPeriod <- function(period, nP)
        layer.apply(seq(1, nP - 12, by = 12), openRatio, period)
    datH = openPeriod(historicID, nHist)
    datF = lapply(futuresID, openPeriod, nFutr)
    mask = !is.na(datH[[1]])

t.test.pSame <- function(v1, v2, ...) {
                if (all(v1 == v2)) return(1)
                t.test(v1, v2, ...)[[3]]
            }
    testFuture <- function(dat, period) {
        dat = addLayer(datH, dat)
        sigFromLast <- function(syr) {
            if (syr < (nlayers(datH) - nsigYrs * 2 - 1)) period = historicID
            tempFile = paste('temp/fireImpact-cells-pdiff', model, variable, period, 
                         syr, nsigYrs , ".nc", sep = '-')
            if (file.exists(tempFile)) return(brick(tempFile))
            print(tempFile)
            index = syr-1 + 1:nsigYrs
            dat1 = FUN(dat[[index]])
            dat2 = FUN(dat[[index + nsigYrs]])
            p0 = 0
            
            sigDif <- function (v1, v2) 
                t.test.pSame(v1, v2, paired = FALSE)
            
            outv = mapply(sigDif, matrix2list(dat1[mask], 1), matrix2list(dat2[mask], 1))
            out = dat1[[1:2]]; out[] = NaN
            out[[1]] = mean(dat1) - mean(dat2)
            out[[2]][mask] = outv
            out = writeRaster(out, file = tempFile, overwrite = TRUE)
            return(out)
        }
        out = layer.apply(1:(nlayers(dat) - nsigYrs*2), sigFromLast)
        browser()
    }
    sigChange = mapply(testFuture, datF, futuresID)     
}

mapply(function(v, FUN) lapply(models, sigTS, v, FUN), variables, FUNs)



