source("libs/filename.noPath.r")
library(ncdf4)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("libs")
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

variables = c("trees", "trees")#, 'cveg', 'csoil') #, '

regionsFile = 'data/GFEDregions.nc'
TCthreshold = 50

gwts = read.csv(gwts)

open2TS <- function(period, experiment, model, variable, mask = NULL, maskID = '') {
    tempFile = paste('temp/fireImpact-GFEDregions-newGWTs', #
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
            
        }
        mapply(forModel, histRatio, rcp)
       
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


outs = lapply(1:14, forRegion)

plot.new()
legend(lty = 1, col = cols, 'top', names(cols), horiz = TRUE, bty = 'n')
legend(lty = 1:4, col = 'black', 'bottom', models, horiz = TRUE, bty = 'n')

dev.off()

sigChangePerRegion <- function(out, name, vt = 2){
    
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
                out = out/mean(out[1])
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
             
            gw = gwts[,c(1, 7, 6, 9, 8, 3, 2, 5, 4)][,1+mod + (rcp-1)*2]
            lines(x, y[, 2], lty = mod, col = "blue")
            lines(x, y[, 1], lty = mod, col = "red")

            findGWT <- function(thresh = 1.5) {
                if (max(gw) < thresh) return(list(thresh, NaN, NaN, NULL, NULL))
                
                threshT = which.min(abs(gw-thresh))
                impact = y[threshT,2]
                yr_noF = x[threshT]
                
                fireTest = diff(y[,1] < impact)
                ## classify no fire
                if (threshT == 1) index = 1:2 else index = (threshT-1):threshT

                testCross <- function(state) {
                    cross = which(fireTest == state)
                    if (length(cross) == 0) {
                        yr_F = threshF = NaN
                        if (diff(y[index,1])*state < 0) fState = 2
                        else {
                            if (any((y[threshT:nrow(y),1] * state) < (y[threshT, 1] * state))) {
                                fState = 2
                            } else {
                                fState = 1
                            }
                        }
                    } else {
                        if (length(cross) > 1) cross = cross[which.min(abs(cross - threshT))]
                        fState = 3
                        yr_F = x[cross]
                        threshF = gw[cross]
                    }
                    return(c(fState, yr_F, threshF))
                }
                
                if (impact < 1) {
                    if (diff(y[index,2]) < 0) {
                        nfState = 1 # reducing
                        c(fState, yr_F, threshF) := testCross(1)
                    } else {
                        nfState = 2 #recovering
                        c(fState, yr_F, threshF) := testCross(-1)                       
                    }
                } else {
                    if (diff(y[index,2]) > 0) {
                        nfState = 3 # increasing
                        c(fState, yr_F, threshF) := testCross(-1)  
                    } else {
                        nfState = 4 #loosing
                        c(fState, yr_F, threshF) := testCross(1)  
                    }
                }

                if (any(thresh == c(1, 1.5, 2))) { 
                    lines(rep(yr_noF, 2), c(0, impact), col = "blue", lty = 2)
                    if (!is.na(yr_F)) {
                        lines(c(yr_noF, yr_F), rep(impact, 2), lty = 2)
                        lines(rep(yr_F, 2), c(0, impact), col = "red", lty = 2)
                    }
                }
                
                return(c(impact, yr_noF, yr_F, thresh, threshF, nfState, fState))
            }
            findGWT(1.5)
            #sapply(c(1, 1.5, 2), findGWT)
            gwi = gw[2:(length(gw)-1)]
            out = sapply(unique(gwi), findGWT)
            #browser()      
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

forRegion <- function(result, regioName, axes = c(), rcp, tFUN) {
    #dev.new() 
    
    plot(c(-0.1, 1), c(-0.1, 1), type = 'n', 
         xaxs = 'i', yaxs = 'i', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
       
        lines(c(-9E9, 9E9), c(-9E9, 9E9))
        if (is.na(result)) {
            x = seq(0.3, 0.9, length.out = length(colsF))

            pnts <- function(xoff = 0, yoff = 0, col = "black") points(x-xoff, x-yoff,  pch = 19, cex = 1.45, col = col)
            pnts(xoff = 0.3)
            pnts(xoff = 0.3, col = colsF)
            pnts(yoff = 0.3)
            pnts(yoff = 0.3, col = colsF)
            text(x = x, y = x-0.3, c("reducing", "recovering", "increasing", "deminishing"), srt = 45, adj = c(0.5, -2), font = c(2, 1, 2, 1))
            text(x = x-0.3, y = x, c("reducing", "recovering", "increasing", "deminishing"), srt = 45, adj = c(0.5, 2), font = c(1, 2, 1, 2))
    
            mtext(side = 1, line = -1, "Earlier", font = 2, adj = 0.67)
            mtext(side = 3, line = -2, "Later", font = 2, adj = 0.33)

            return()
        }
        mtext(side = 3, line = -2, adj = 0.1, regioName)
        
        labels = c(0.2, 0.5, 1, 1.5, 2, 3, 4)
        labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
        at = tFUN(labels)
        labels[length(labels)] = '3.5+'
        axisFUN <- function(side) axis(side, at = at, labels = labels)
        lapply(axes, axisFUN)
        labels = c(0.2, 0.5, 1, 1.5, 2, 3, 4)
        labels = c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
        at = tFUN(labels)
        labels[length(labels)] = '3.5+'
        axisFUN <- function(side) axis(side, at = at, labels = labels)
        lapply(axes, axisFUN)
        xg = seq(-10, 10, 0.01)
        for (dy in seq(-4, 4, 0.5)) lines(tFUN(xg), tFUN(xg-dy), lty = 4)
        lines(tFUN(c(3.5, 3.5)), c(-9E9, 9E9));lines(c(-9E9, 9E9), tFUN(c(3.5, 3.5)))
        #forRCP <- function(rcp, pch) {
        
        forModel <- function(mod) {
            
            res = result[[rcp]][[mod]]
            res[4:5,] = tFUN(res[4:5,])
            
            points(res[4,], res[5,], pch = 19, cex = 1.45)
            for (cex in seq(1.2, 0.1, -0.1))
                points(res[4, ], res[5,], col = colsF[res[6,]], pch = 19, cex = cex)
            test = is.na(res[5,])
            off = res[6,test]; on = res[7, test];
            
            pchT = ((res[6,test] == 1 | res[6,test] == 4) & 
                    (res[7,test] == 1 | res[7, test] == 4)) | 
                   ((res[6, test] == 2 | res[6, test] == 3) & 
                    (res[7, test] == 2 | res[7, test] == 3))
            yexceed = seq(0.9, 1, length.out = sum(test))
            yexceed = sample(yexceed, sum(test), replace = FALSE)
            pch = c(4, 19)[pchT+1]
            points(res[4,test], yexceed, pch = pch, cex = 1.45)
            points(res[4, test], yexceed,pch = pch, cex = 1.2, col = colsF[res[6,test]])
            return()
        }
        
        mapply(forModel, 1:4)#, c("#1b9e77", "#d95f02", "#b2df8a", "#a6cee3"))
        
}

tFUN <- function(x) {
    test = x > 3.5
    x[test] = 4
    x/4# {y = (1-exp(-abs(x)))^2; y[x<0] = -y[x<0]; y}
}

colsF = c("#8c510a", "#c7eae5", "#01665e", "#f6e8c3")
colsT = c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7',
          '#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')
plot4RCP <- function(rcp){
    png(paste0("figs/GWTs_new", rcp, ".png", sep = '-'), 
        height = 12, width = 12, res = 300, units = 'in')
        lmat = t(matrix(c(1:15, 0), ncol = 4))
        layout(lmat)
        par(oma = c(3, 3.5, 1, 2), mar = rep(0.5, 4))
        axis3 = c(rep(T, 4), rep(F, 10))
        axes = cbind(rev(axis3), rep(c(T, F, F, F), length.out = 14), 
                    axis3, rep(c(F, F, F, T), length.out = 14))
    
        mapply(forRegion, results, regionNames, axes = apply(axes, 1, which), rcp, 
               MoreArgs = list(tFUN = tFUN))
        mtext.units("Without fire (~DEG~C)", side = 1, outer = TRUE, line = 2)
        mtext.units("With fire (~DEG~C)", side = 2, outer = TRUE, line = 2)
        par(mar = c(0.5, 2, 2, 0.5))
        forRegion(NaN, NaN, NaN, NaN, tFUN)
    dev.off()
}
mapply(plot4RCP, 1:2)


parisSumm <- function(gwt = 1.5) {
    forRegion <- function(res) {
        forRCP <- function(xs) {
            forMod <- function(x) {
                          
                index = which.min(abs(gwt - unlist(x[4,])))
                out = outi = x[4:7, index]
                out[1:2] = round(out[1:2], 2)
                #if (is.na(out[2])) browser()
                out[3] = c("red", "rec", "inc", "dem")[out[3]]
                out[2] = c("WD", "NR", out[2])[outi[4]]
                return(out[1:3])
            }
            out = sapply(xs, forMod)
        }
        out = lapply(res, forRCP)
        out = do.call(rbind, out)
        
        rownames(out) = paste(rep(c('RCP2.6', 'RCP6.0'), each = 3), c('NoFire', 'WithFire', 'Pathway'))
        return(out)
    }
    out = lapply(results, forRegion)
    out = do.call(rbind, out)
    
    rownames(out) = paste(rep(regionNames, each = 6), rownames(out), sep = '-')
    return(out)
}

out = lapply(c(1, 1.5, 2), parisSumm)
out = do.call(cbind, out)

colnames(out) = paste('GWT', rep(c(1, 1.5, 2), each = 4), models, sep = '-')


write.csv(out, "outputs/GWL_equiv.csv")
browser()




browser()
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
