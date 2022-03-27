source("libs/filename.noPath.r")
library(ncdf4)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
source("libs/process_jules_file.r")
source("libs/writeRaster.Standard.r")
options(error=recover)
graphics.off()

dir = '../ConFIRE_ISIMIP/INFERNOonOff/Global/'
gwts = 'data/Tas_vs_year-ISIMIP.csv'
historicID = "historic"
futuresID = c("RCP2.6", "RCP6.0")

experiments = c("_on", "_off")
models = c("GFDL-ESM2M", "HADGEM2-ES", "MIROC5", "IPSL-CM5A-LR")

variables = c("trees", 'cveg', 'csoil') #, '

regionsFile = 'data/GFEDregions.nc'
TCthreshold = 50

gwts = read.csv(gwts)
gwts = gwts[, 1+c(0, 6, 5, 8, 7, 2, 1, 4, 3)]
missing = sapply( gwts[1, -1], '*', ((1871:(gwts[1,1]-1))-1871)/(gwts[1,1]-1871))
missing = cbind(1871:(gwts[1,1]-1), missing)
gwts = rbind(missing, as.matrix(gwts))

open2TS <- function(period, experiment, model, variable, mask = NULL, maskID = '') {
    tempFile = paste('temp/fireImpact-GFEDregions', 
                     period, experiment, model, variable, maskID, 
                     'time series.Rd', sep = '-')
    print(tempFile)
    if (file.exists(tempFile)) {
        load(tempFile)
        return(out)
    }
    dat = brick(paste0(dir, '/', period, experiment, '/', model, '/', variable, '.nc'))
    
    mask = round(raster::resample(mask, dat[[1]]))
    maskA = raster::area(mask, na.rm  = TRUE)[mask == 1]
    
    sindex = seq(1, nlayers(dat), 50)
    
    maskSum <- function(si) {
        print(si)
        r = dat[[si:min(si+49, nlayers(dat))]]
        apply(r[mask == 1,], 2, function(vs) sum(vs*maskA, na.rm = TRUE))
    }
    out = lapply(sindex, maskSum)
    out = do.call(c, out)
    
    if (period == "historic" && variable == "trees" && experiment == "_on" && grepl( 'All', maskID)) {
        TC = mean(dat[[1:480]])
        tcMask = addLayer(TC>0.4 & mask, TC<0.4 & mask)
        out = list(out, tcMask)
    } else out = list(out, NULL)
    save(out, file = tempFile)
    return(out)
}

forVariable <- function(variable, masks, maskID, ...) {
    forModel <- function(model, mask) {
        forExperiment <- function(experiment) {
            hist = open2TS(historicID, experiment, model, variable, mask = mask, maskID = maskID, ...)
            futr = lapply(futuresID, open2TS, experiment, model, variable, mask = mask,maskID = maskID, ...)
            futr = lapply(futr, function(i) i[[1]])
            tcMask = hist[[2]]
            hist = hist[[1]]
            return(list(hist, futr, tcMask))
        }
        out = lapply(experiments, forExperiment)
        
        list(sapply(out, function(i) i[1:2]), out[[1]][[3]])
    }
    dat = mapply(forModel, models, masks, SIMPLIFY = FALSE)
    plotVar(dat, variable, maskID) 
    
    return(dat)
}

plotVar <- function(dat, nme, maskID) {
    dat = lapply(dat, function(i) i[[1]])
    
    histOn = lapply(dat, function(i) i[[1,1]])
    histOff = lapply(dat, function(i) i[[1,2]])

    openFut <- function(id) {
        out = lapply(dat, function(i) i[[2,id]])
        lapply(1:2, function(i) lapply(out, function(j) j[[i]]))
    }
    futrOn = openFut(1)
    futrOff = openFut(2)#lapply(dat, function(i) i[[2,2]])
    
    xhist = seq(1861, by= 1/12, length.out = length(histOn[[1]]))
    xfutr = seq(2005, by= 1/12, length.out = length(futrOn[[1]][[1]]))
    xrange = range(c(xhist, xfutr))
    yrange = range(c(histOn, histOff, unlist(futrOn), unlist(futrOff)))
    plot(xrange, yrange, type = 'n', xlab = '', ylab = '')

    labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)
    at = yrange[1] + diff(yrange) *labels
    axis(4, at = at, labels = labels)
    plotLines <- function(xs, ys,col = 'black')
        mapply(lines, ys, lty = 1:4, MoreArgs = list(x = xs, col = col))

    
    plotLines(xhist, histOn, cols['histOn'])
    plotLines(xfutr, futrOn[[1]], cols['futrOn26'])
    plotLines(xfutr, futrOn[[2]], cols['futrOn60'])
   
    plotLines(xhist, histOff, cols['histOn'])
    plotLines(xfutr, futrOff[[1]], cols['futrOff26'])
    plotLines(xfutr, futrOff[[2]], cols['futrOff60'])

    divideExpand <- function(i,j){
        out = i/j
        yrange[1] + diff(yrange) * out
    }

    histRatio = mapply(divideExpand, histOn, histOff, SIMPLIFY = FALSE)
    #mapply(lines, histRatio, lty = 1:4, MoreArgs = list(x = xhist))

    RCP26Ratio = mapply(divideExpand, futrOn[[1]], futrOff[[1]], SIMPLIFY = FALSE)
    RCP60Ratio = mapply(divideExpand, futrOn[[2]], futrOff[[2]], SIMPLIFY = FALSE)
    
    #  mapply(lines, histRatio, lty = 1:4, MoreArgs = list(x = xhist))
    plotLines(xhist, histRatio, cols['histRatio'])
    plotLines(xfutr, RCP26Ratio, cols['RCP26Ratio'])
    plotLines(xfutr, RCP60Ratio, cols['RCP60Ratio'])
    
    forRCP <- function(rcp) {
        forModel <- function(h, f) {
            full = c(h, f)
            forPeriod <- function(end) {
                p1 = full[(end-479):(end-240)]
                p2 = full[(end-239):end]
                t.test(p1, p2)[[3]]
            }
            sig = sapply(480:length(full), forPeriod)
            browser()
        }
        mapply(forModel, histRatio, rcp)
        browser()
    }
    #forRCP(RCP26Ratio)  
    mtext(side = 3, paste(nme, '-', maskID))
    return(list(histRatio, RCP26Ratio, RCP60Ratio))
}

regions = raster(regionsFile)
forRegion <- function(id) {
    mask = regions == id
    dats = lapply(variables, forVariable,   masks = rep(list(mask), length(models)),
                  maskID = paste0(id, '-All'))
   

    forMask <- function(mno, msname = '') {
        masks = lapply(dats[[1]], function(i) i[[2]][[mno]])
        dats = lapply(variables, forVariable,   masks = masks,
                    maskID = paste0(id, msname))
    }
    datsF = forMask(1, '-forest')
    datsN = forMask(2, '-nonforest')
    return(list(dats, datsF, datsN))
}

png("figs/TS-stuff.png", height = 7, width = 9, res = 300, units = 'in')
layout(rbind(matrix(1:(length(variables)*3), nrow = 3), 1+length(variables)*3))
par( mar = c(2, 2.5, 2, 2.5), oma = c(3, 0, 1, 0))
cols = c(histOn = '#1b9e77', futrOn26 = '#7570b3', futrOn60 = '#d95f02',
             histOff = '#66c2a5', futrOff26 = '#8da0cb', futrOff60 = '#fc8d62',
             histRatio = 'black', RCP26Ratio = 'blue', RCP60Ratio = 'red')


#outs = lapply(1:14, forRegion)

plot.new()
legend(lty = 1, col = cols, 'top', names(cols), horiz = TRUE, bty = 'n')
legend(lty = 1:4, col = 'black', 'bottom', models, horiz = TRUE, bty = 'n')

dev.off()


sigChangePerRegion <- function(out, name, vt = 1){
    out = out[[vt]][[2]]
    
    forModel <- function(mout) {
        mout = mout[[1]]
        
        forRCP <- function(rcp) {
            dat = apply(mout, 2, function(i) c(i[[1]], i[[2]][[rcp]]))
            syrs = seq(1, nrow(dat) - 11, 12)
            annualSum <- function(syr, x, dt = 11) 
                sum(x[syr:(syr+dt)])
            
            annualSums <- function(...) sapply(syrs, annualSum, ...)

            annualSumsRatio <- function(..., nyrs = 20) {
                out = annualSums(...)
                out = sapply(1:(length(out)-20), annualSum, out, 20)/20
                out = out/mean(out[131])
            }   
            dat = apply(dat, 2, annualSumsRatio)
        }
        lapply(1:2, forRCP)
    }
   
    ratios = lapply(out, forModel)
    
    x = 1:nrow(ratios[[1]][[1]]) + 1870
    
    forRCP <- function(rcp) {
        rcpR = lapply(ratios, function(i) i[[rcp]])
        
        plot(range(x), range(unlist(rcpR)),xaxt = 'n', 
             type = 'n', xlab = 'years', ylab = 'ratio', yaxt = 'n')
        if (name == tail(regionNames, 1)) axis(1)
        if (name ==regionNames[1]) axis(3)
        if (rcp == 1) axis(2)
        if (rcp == 2) axis(4)
        grid()
        mtext(name, side = 2, adj = 0.1, line = -2)    
        plotMod <- function(mod) {
            y = rcpR[[mod]]
            gw = gwts[,1+mod + (rcp-1)*2]
            lines(x, y[, 2], lty = mod, col = "blue")
            lines(x, y[, 1], lty = mod, col = "red")

            findGWT <- function(thresh = 1.5) {
                if (max(gw) < thresh) return(list(thresh, NaN, NULL, NULL))
                
                index = which.min(abs(gw-thresh))

                yr = x[index]
                val_noFire = y[index, 2]
                val_Fire = y[index, 1]

                xeq <- function(a, b) if (val_noFire < 1) return(a <= b) else return(a >= b)
                lines(c(yr, yr), c(0, val_noFire), lty = mod)

                find4Experiment <- function(id) {
                    index = range(which(xeq(y[,id], val_noFire)))
                    yr = x[index]
                    gwts = gw[index] 
                    return(list(yr = yr, gwts = gwts))
                }
                noFire = find4Experiment(2)
                Fire = find4Experiment(1)
                

                lines(noFire[['yr']], c(val_noFire, val_noFire), lty = mod)
                lines(noFire[['yr']], c(val_noFire, val_noFire) + diff(par("usr")[3:4])*0.002,
                      lty = mod)
            
                lines(Fire[['yr']], c(val_noFire, val_noFire), lty = mod)
                
                            
                return(list(thresh,yr,  noFire, Fire))
            }
            out = cbind(findGWT(1.5), findGWT(2.0))
            return(out)
        }
        out = mapply(plotMod, 1:4, SIMPLIFY = FALSE)

        return(out)
    }
    outs = lapply(1:2, forRCP)
    return(outs)
}
regionNames = c("BONA", "TENA", "CEAM", "NHSA", "SHSA", "EURO", "MIDE", "NHAF", "SHAF", "BOAS",
                "CEAS", "SEAS", "EQAS", "AUST") 
png("figs/fireGWTs_cal.png", height = 42, width = 10, res = 300, units = 'in')
    par(mfrow = c(length(outs),2), mar = c(0.25, 2.5, 0.25, 0.25), oma = c(3, 3, 0, 0))
    results = mapply(sigChangePerRegion, outs, regionNames,
                     MoreArgs = list(2), SIMPLIFY = FALSE)
    mtext(side = 2, outer = TRUE, "With:without fire")
dev.off()

forRegion <- function(result, regioName, axes = c()) {
   
    tFUN <- function(x) {y = (1-exp(-abs(x)))^2; y[x<0] = -y[x<0]; y}
    plot(c(0, 1), c(0, 1), type = 'n', 
         xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    mtext(side = 3, line = -2, adj = 0.1, regioName)

    labels = c(0.2, 0.5, 1, 1.5, 2, 3, 4)
    axisFUN <- function(side) axis(side, at = tFUN(labels), labels = labels)
    lapply(axes, axisFUN)
    
    lines(c(-9E9, 9E9), c(-9E9, 9E9))

    xg = seq(-10, 10, 0.01)
    for (dy in seq(-4, 4, 0.5)) lines(tFUN(xg), tFUN(xg-dy), lty = 4)

    forRCP <- function(rcp, pch) {
        forTemp <- function(tp, col) {
        
            temps = lapply(result[[rcp]], function(i) sapply(i[, tp][3:4], function(j) j[[2]]))
            errorBox <- function(xy)  {
                if (all(sapply(xy, is.null))) return()
                cp = apply(xy, 2, mean)
                #arrows(xy[1,1], cp[2], xy[2,1], cp[2], angle = 90, code = 3, col = col)
                #arrows(cp[1], xy[1,2], cp[1], xy[2,2], angle = 90, code = 3, col = col)
                #points(cp[1], cp[2], pch  = pch, col = col)
                if (any(is.na(xy))){
                    print("yay")
                    xy[is.na(xy)] = 1000
                }
                points(tFUN(xy[1,1]), tFUN(xy[1,2]), pch = pch, col = col)
            }
            lapply(temps, errorBox)
            
        }
        mapply(forTemp, 1:2, c("blue", "red"))
    }
    mapply(forRCP, 1:2, c(4, 20))
}
png("figs/GWTs_new.png", height = 12, width = 12, res = 300, units = 'in')
    par(mfrow = c(4, 4), oma = c(2, 2, 1, 2), mar = rep(0.5, 4))
    axis3 = c(rep(T, 4), rep(F, 10))
    axes = cbind(rev(axis3), rep(c(T, F, F, F), length.out = 14), 
                 axis3, rep(c(F, F, F, T), length.out = 14))
    
    mapply(forRegion, results, regionNames, axes = apply(axes, 1, which))
dev.off()

summerize <- function(dats) {
    forRCP <- function(rcp) {
        out = apply(sapply(dats, function(i) i[,rcp]), 1, range)
        c(paste(out[,1], collapse = '-'), paste(round(out[,2], 2), collapse = '-'))
    }
    
    out = sapply(1:2, forRCP)
    colnames(out) = c("RCP2.6", "RCP6.0")
    rownames(out) = c("year", "GWT")
    out
}
browser()
out = sapply(results, summerize)
colnames(out) = regionNames
write.csv(out, file = "yay2.csv")
