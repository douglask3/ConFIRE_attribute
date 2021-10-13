library(raster)
library(rasterExtras)
library(gitBasedProjects)
library(ncdf4)
library(mapproj)
sourceAllLibs("libs/")
graphics.off()

cols = rev(c("#01665e", "#8c510a", 'red'))
fname = "figs/burnt2021/"
grabCache = FALSE

cveg = read.csv("data/Burton2021_cveg.csv", stringsAsFactors = FALSE)[-1,]
cveg = spline(cveg[,1], cveg[,2], n = 1000)

fColGen <- function(col) {
    xfCols = c("white", rep(col, 10))
    xfCols =  make_col_vector(xfCols, limits = 1:256)
    hx = c(0:9, 'A', 'B', 'C', 'D', 'E', 'F')
    hx = rev(paste0(rep(hx, each = length(hx)), hx))
    hx = hx[round(256*seq(0, 1, length.out = 256)^0.3)]
    xfCols = paste0(xfCols, hx)
}

xfCols = fColGen('blue')
yfCols = fColGen('red')

tas = read.csv("data/Burton2021_tas.csv", stringsAsFactors = FALSE)[-1,]  
tas = spline(tas[,1], tas[,2], n = 1000) 

tasYrs <- function(s, e, x) {
    test = x[[1]] > s & x[[1]] < e
    mean(x[[2]][test])
}
yrs = 1890:2100

mnVeg = mapply(tasYrs, yrs-30, yrs, MoreArgs = list(cveg))
cveg[[2]] = cveg[[2]] - max(mnVeg)
limits = seq(min(cveg[[2]]), max(cveg[[2]]), length.out = 20)
cols =  make_col_vector(cols, limits = limits)


mnTas = mapply(tasYrs, yrs-30, yrs, MoreArgs = list(tas))
GWT = c(1.5, 2, 3, 4, 5)
GWT = cbind(yrs[sapply(GWT, function(i) which(mnTas > i)[1])], paste0(GWT, '~DEG~C'))
YrHts = rbind(c(2021, 'Today'), GWT)

# ourworldindata,
CO2_2_C = 12/((16*2) + 12) 
CO2e_2019 = 36.44 * CO2_2_C
UK_acc2019 = 77.84 * CO2_2_C
CH_acc2019 = 219.99 * CO2_2_C
US_acc2019 = 410.25 * CO2_2_C/2
             
cveg2day =  mnVeg[yrs == 2021]
dcveg2day = cveg[[2]][which.min(abs(cveg[[1]] - 2021))]
whichCyr <- function(x)  {
    x = dcveg2day - x
    test = which(cveg[[2]] < x)[1]
    #test = which((cveg2day - mnVeg) > x)[1]
    c(cveg[[1]][test], x)
}
CrHts = rbind(c(cveg[[1]][which(abs(cveg[[2]]) > cveg2day/2)[1]], -cveg2day/2, 'Half Amazon carbon', '1/2 Amazon'),
             c(whichCyr(CO2e_2019), '2019 global ~CO2~e emission', '2019 global emissions'),
             c(whichCyr(UK_acc2019), 'UK cumulative  ~CO2~e emissions 1750-2019', 'UK'),
             c(whichCyr(CH_acc2019), 'China cumulative  ~CO2~e emissions', 'China'),
             c(whichCyr(US_acc2019), 'Half USA cumulative  ~CO2~e emissions', 'USA/2'))
#browser()
plotFun <- function(stYr, endYr, flash, plotN) {
    xlim = c(stYr, endYr)
    plotN = paste(c(rep('0', 6-nchar(plotN)), plotN), collapse = '')
    fname = paste0(fname, '/plot-', plotN, '.png')
    if (grabCache && file.exists(fname)) return()

    
    test = cveg[[1]] >= xlim[1] & cveg[[1]] <= xlim[2]
    x = cveg[[1]][test]; y = cveg[[2]][test]
    col = cols[cut_results(y, limits)]   

    png(fname, height = 5, width = 5, res = 300, units = 'in')
    test = cveg[[1]] >= xlim[1] & cveg[[1]] <= xlim[2]
    x = cveg[[1]][test]; y = cveg[[2]][test]
    col = cols[cut_results(y, limits)]   
    
    plot(x, y, type = 'n', axes = FALSE, xlab = '', ylab = '')
    polygon(9E9 * c(-1, 1, 1, -1), 9E9 * c(-1, -1, 1, 1), col = 'black', xpd = NA)
    

    addFlash <- function(Hts, s, e, cols) {
        Ht = which(Hts[,1] > s & Hts[,1] < e)
        if (length(Ht) > 0) {
            dist = round(256*(min(e - as.numeric(Hts[Ht,1]))/ (e - s)))
            col = cols[dist]
            polygon(9E9 * c(-1, 1, 1, -1), 9E9 * c(-1, -1, 1, 1), col = col)
        } 
    }
    if (flash) {
        addFlash(YrHts, stYr, endYr,xfCols) 
        addFlash(CrHts, stYr, endYr,yfCols) 
    }
    addSeg <- function(i) 
        lines(x[(i-1):i], y[(i-1):i], col = col[i], lwd = 5)

    lapply(1:length(y), addSeg)
    #xHt = which(YrHts[,1] > stYr & YrHts[,1] < endYr)
    #if (length(xHt) > 0) {
    #    dist = round(256*(min(endYr - as.numeric(YrHts[xHt,1]))/ (endYr - stYr)))
    #    flCol = xfCols[dist]
    #    polygon(9E9 * c(-1, 1, 1, -1), 9E9 * c(-1, -1, 1, 1), col = flCol)
    #}

    addHt <- function(x) {
        lines(rep(as.numeric(x[1]), 2), c(-9E9, 9E9), lwd = 2, col = 'white')
        text.units(x[2], x = as.numeric(x[1]), y = max(y), adj = -0.5, srt = 90, 
                   col = 'white', xpd = NA)
    }
    apply(YrHts, 1, addHt)
    x0 = x
    addHt <- function(x) {
        lines(c(-9E9, 9E9), rep(as.numeric(x[2]), 2), lwd = 2, col = 'white')
        if (flash) str = x[3] else str = x[4]
        if (nchar(str) > 6) adj = c(1, 2) else adj = c(-0.1, 0)
        
        text.units(str, y = as.numeric(x[2]), x = endYr, adj = adj, col = 'white', xpd = NA)
    }
    apply(CrHts, 1, addHt)
    
    yrTick = unique(round(seq(xlim[1], xlim[2], length.out = 5)))
    axis(2, lwd = 3, font = 2, col = 'white', col.axis = 'white')
    axis(1, at = yrTick, lwd = 3, font = 2, col = 'white', col.axis = 'white')
    
    dev.off()
}

makeDir(fname)

start = 1861:2070
end = start + 30

start = c(rep(1861, 30), start)
end = c(1862:1891, end)

start = 1890:2070
start = c(start, rep(2070, 5))
flash = rep(T, length(start))
start = c(rep(start[1], 20), start)
end = start + 30

start = c(start, seq(2070, 1860, length.out = 20))
end = c(end, rep(2100, 20))

start = c(start, rep(1860, 50))
end = c(end, rep(2100, 50))
flash = c(flash, rep(F, length(start) - length(flash)))

mapply(plotFun, start, end, F, 1:length(start))
commd = paste0("convert -delay 10 ", fname, '/*.png ', fname, '/out.gif')
system(commd)

