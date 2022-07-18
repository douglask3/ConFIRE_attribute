graphics.off()

dat = read.csv("data/Tas_vs_year.xlsx - NBP.csv", header = FALSE, stringsAsFactors = FALSE)

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
        idM = which(dat[,2] == model)
        id = idR[idR > idM][1]
        apply(dat[id+(2:1), -1], 2, as.numeric)
    }

    nbps = lapply(models, modelNBP)
    dnps = lapply(nbps, function(i) i[2,] - i[1,])
    
    testKicker <- function(x, temp, syrs = 2000) {
        syr = which(years == syrs)
        x0 = x[syr+ 1:20]
        doTest <- function(i) 
            wilcox.test(x0, dnps[[1]][i + 1:20], paired = FALSE)[[3]]

        pvs = sapply(syr:(length(x) -20), doTest)
        sig = which(pvs < 0.1)[1]
        c(sig + syrs + 10, temp[sig + syr])
    }
    fireSig = mapply(testKicker, dnps, temps, SIMPLIFY = FALSE)
    
    findColsLims <- function(dats) {
        limits = sort(c(0, quantile(dats, seq(0.1, 0.9, 0.1))))
        limits = unique(signif(limits, 1))
    
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
   
    
    forModel <- function(model, mi, x, mark) {
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
                plotLines(mark[1], y[1], 'black', lwd = 1.5, lty = 2)
                text(x = mark[1], y = y[1], adj = c(0.5, -0.2), round(mark[2], 2), srt = 90)
            }
        }
        forExp(1, cols, limits)        
        forExp(2, cols, limits)      
        forExp(3, dcols, dlimits)
    }
    mapply(forModel, models, 1:length(models), xs, fireSig)

    legendColBar <- function(xx, yy, dxl, cols, limits, switch = FALSE,
                             extend_max = TRUE, extend_min = TRUE) {

        ys = seq(yy[1], yy[2], length.out = length(cols) +1)
        addBox <- function(y1, y2, yi, col) {
            print(yi)
            polyFun <- function(x, y) polygon(x, y, col = col, lwd = 2, xpd = TRUE)
            if (yi == 1 && extend_min) 
                polyFun(c(xx, mean(xx), xx[1]), c(y2, y2, y1, y2))
            else if (yi == length(cols) && extend_max) 
                polyFun(c(xx, mean(xx), xx[1]), c(y1, y1, y2, y1))
            else polyFun(c(xx, rev(xx), xx[1]), c(y2, y2, y1, y1, y2))
                 
        }
        id = 1:length(cols)
        if (switch) {
            id = rev(id)
            cols = rev(cols)
            limits = rev(limits)
        }
        mapply(addBox, ys[-1], head(ys, -1), id, cols)
        text(xx[2] + diff(xx) *0.3, head(ys[-1], -1), limits, xpd = NA, adj = 0)
    }
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
