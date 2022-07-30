###########
## setup ##
###########
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
#library(rasterExtras)
#library(gitBasedProjects)
library(ncdf4)
sourceAllLibs("libs/")
graphics.off()

dat = read.csv("data/Tas_vs_year.xlsx - NBP.csv", header = FALSE, stringsAsFactors = FALSE)#[,-1]

models = c("HadGEM2-ES" = "H", "GFDL" = "G", "IPSL" = "I", "MIROC" = "M")

regions = c("Global", "BONA", "TENA", "CEAM", "NHSA", "SHSA", "EURO", "MIDE", 
                      "NHAF", "SHAF", "BOAS", "CEAS", "SEAS", "EQAS", "AUST")



modTemp <- function(model) {
    id = which(dat[,2] == model) + 1
    if (model == 'H') id = id + 1
    out =  as.numeric(dat[id,-1])
    out[1] = 0
    out
}
temps = lapply(models, modTemp)

cols1 = c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5')
cols2 = c('#c7eae5','#80cdc1','#35978f','#01665e','#003c30')

years = as.numeric(dat[2,-1])
years[1] = years[2] - 1
forRegions <- function(region, Rid, yearsTest = TRUE) {
    if (region == "Global") idR = which(dat[,1] == "GlobTemp")
    else idR = which(dat[,2] == region)
    
    if (yearsTest) xs = lapply(temps, function(i) years)
    else xs = temps#lapply(temps, function(i)  i[tyrID,-1])
    xrange = range(unlist(xs))
    
    plot(xrange, c(0, 14), axes = FALSE, 
         ylim = c(13, 0), xlab = '', ylab = '', type = 'n')
    
    #axis(1)
    modelNBP <- function(model) {
        idM = which(as.character(dat[,2]) == model)
        id = idR[idR > idM][1]
        
        apply(dat[id+(2:1), -1], 2, as.numeric)
    }

    nbps = lapply(models, modelNBP)
    dnps = lapply(nbps, function(i) i[2,] - i[1,])
    
    testKicker <- function(x, temp, syrs = 2000, offset = 10) {
        syr = which(years == syrs)
        x0 = x[syr+ 1:20]
        doTest <- function(i) 
            wilcox.test(x0, dnps[[1]][i + 1:20], paired = FALSE)[[3]]

        pvs = sapply(syr:(length(x) -20), doTest)
        sig = which(pvs < 0.1)[1] + offset
        c(sig + syrs + 10, temp[sig + syr])
    }
    fireSig = mapply(testKicker, dnps, temps, SIMPLIFY = FALSE)
    
    ddnps = lapply(dnps, function(x) sapply(20:length(x), function(i) mean(x[(i-19):i])))
    fireSig_an = mapply(testKicker, ddnps, temps, SIMPLIFY = FALSE)
    
    findColsLims <- function(dats) {
        limits = sort(c(0, quantile(dats, seq(0.1, 0.9, 0.1))))
        limits = unique(signif(limits, 1))
        limits[limits == -0.09] = -0.1
        limits = limits[limits != -0.0003]
 
        colsA = colsB = c()
        if (any(limits < 0))
            colsA = make_col_vector(cols1, ncols = 1+sum(limits <= 0))
        if (any(limits > 0))
            colsB = make_col_vector(cols2, ncols = 1+sum(limits >= 0))
        cols = c(colsA, colsB)
        if ((length(cols) -2) == length(limits)) cols = c(colsA, colsB[-1])
        if ((length(cols) -3) == length(limits)) cols = c(head(colsA, -1), colsB[-1])
        if ((length(cols) -1) != length(limits)) browser()
        return(list(limits, cols))
    }
    c(limits, cols) := findColsLims(unlist(nbps))
    c(dlimits, dcols) := findColsLims(unlist(dnps))
   
    forModel <- function(model, mi, x, mark1, mark2) {
        idM = which(dat[,2] == model)
        id = idR[idR > idM][1]
        nbp = dat[id+(2:1), ]
        forExp <- function(i, cols, limits) {
            y = mi + (i-1) * (length(models)+0.5) 
            y = rep(y, length(x))
            if (i == 3) col = as.numeric(nbp[2,-1]) - as.numeric(nbp[1,-1])
            else col = as.numeric(nbp[i,-1])
            col = cut_results(col, limits) 
            plotLines <- function(x, y, col, lwd = 3, ...) 
                lines(c(x, x), y + c(-0.5, 0.5), col = col, lend = 1, lwd = lwd, ...)
            cols.pt = cols[col]  
            mapply(plotLines, x, y, cols.pt)
            if (i == 3) {
                addMark <- function(mark, lty = 2, col = "black", adj = -0.2, srt = 270) {
                    plotLines(mark[1], y[1], col, lwd = 1.5, lty = lty)
                    text(x = mark[1], y = y[1], col = col,
                        adj = c(0.5, adj), round(mark[2], 2), srt = srt)
                }
                addMark(mark1, 2); addMark(mark2, 3, "white", adj = -0.2, srt = 90)
            }
        }
        forExp(1, cols, limits)        
        forExp(2, cols, limits)      
        forExp(3, dcols, dlimits)
    }
    mapply(forModel, models, 1:length(models), xs, fireSig, fireSig_an)
    
    xx = xrange[2] + diff(xrange) * c(0.017,0.05)
    legendColBar(xx, c(0.5, 8.5), 10, cols, limits, TRUE)
    legendColBar(xx, c(9, 14), 10, dcols, dlimits, TRUE)
    mtext(side = 3, adj = 0.1, region, line = -1.7)
    if (Rid %% 4 == 1) {
        textSt <- function(y, txt) 
            text(x = xrange[1] - diff(xrange)*0.04, y = y , txt, srt = 90, cex = 1.3, xpd = NA)
        
        textSt(02.5 , 'Without Fire')
        textSt(07.0 , 'With Fire')
        textSt(11.5, 'Difference')
    }
    axis(1)
   # if (Rid > (length(regions) -4)) axis(1)
}

png("figs/NBPstripes.png", res = 300, units = 'in', height = 14, width = 12)
par(mfrow = c(4, 4), mar = c(1.5, 1.5, 0, 1.5), oma = c(1.5, 1, 1, 1.5))
mapply(forRegions, regions, 1:length(regions), TRUE)
dev.off()
