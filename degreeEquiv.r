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
futuresID = "RCP6.0"

experiments = c("_off", "_on")
models = c("GFDL-ESM2M", "HADGEM2-ES", "MIROC5", "IPSL-CM5A-LR")

variable = "trees"#, 'cveg', 'csoil') #, '

tfile0 = 'temp/degreeEquiv'

openDat <- function(period, experiment, model, variable) 
    brick(paste0(dir, period, experiment, '/', model, '/', variable, '.nc'))


aa <- function(i, dat, ..., ncount = 11) {
    print(i)
    tfile = paste0(c(tfile0, ...,  ncount, i, '.nc'), collapse = '-')
    if (file.exists(tfile)) return(raster(tfile))
    print(tfile)
    out = mean(dat[[(i-ncount):i]])
    writeRaster(out, file = tfile, overwrite = TRUE)
}

openAllP <- function(...) {
    tfile = paste0(c(tfile0, 'running21all', futuresID,  ..., '.nc'), collapse = '-')
    if (file.exists(tfile)) return(brick(tfile))
    hist = openDat(historicID, ...)
    futr = openDat(futuresID, ...)
    dat = addLayer(hist, futr)
    datYr = layer.apply(seq(12, nlayers(dat), by = 12), aa, dat, ...)
    dat = layer.apply(21:nlayers(datYr), aa, datYr, 'running21', ..., ncount = 20)
    dat = dat/dat[[1]]
    dat = writeRaster(dat, file = tfile, overwrite = TRUE)
}

openMod <- function(...) 
    dat = lapply(experiments, openAllP, ...)


dats = lapply(models, openMod, variable = variable)
#browser()

analyseCell <- function(i, dat, gwt, Temp = 1.5, ID) {
    ts1 = as.vector(dat[[1]][i]); ts2 = as.vector(dat[[2]][i])
    if (all(ts1 == 0) && all(ts2 == 0) ) return(rep(NaN, 4))

    tgrad  <- function(id, ts) diff(ts[(id-1):id]) > 0
    grad = tgrad(ID, ts1)
    impact = ts1[ID]
    tt2 = (ts2 - impact) > 0
    f_ID = 1+which(tt2[-1] != head(tt2, -1))
    #f_ID = which.min(abs(ts2 - impact))
    if (length(f_ID) == 0) f_ID = length(ts2)

    f_grad = sapply(f_ID, tgrad, ts2)
    
    if (length(f_ID) > 1) {
        gtest = which(f_grad == grad)
        if (length(gtest) > 0) {
            f_ID = f_ID[gtest]
            f_grad = f_grad[gtest]
        } 
        if (length(f_ID)>0) {
            gtest = which.min(abs(f_ID-ID))
            f_grad = f_grad[gtest]
            f_ID = f_ID[gtest]
        }
    }
    if (f_ID == 1) browser()
    return(c(impact, grad, f_ID, f_grad))
}

gwts = read.csv(gwts)[-1,]


analyseModel <- function(dat, model, Temp = 1.5, switch = FALSE) {
    if (switch) tfile = paste('outputs/degreeEquiv-switch', model, Temp, '.nc', sep = '-')
    else tfile = paste('outputs/degreeEquiv', model, Temp, '.nc', sep = '-')
    if (file.exists(tfile)) return(brick(tfile))
    mask = !is.na(dat[[1]][[1]])
    gwt = gwts[, grepl(substr(futuresID, 4, 6), names(gwts)) & 
                 grepl(substr(model, 1, 4), names(gwts), ignore.case = TRUE)]
    if (switch) dat = dat[2:1]
    index = which(mask[])#[1:42000]
    ncells = length(index)
    ntest = 1000
    analyseGroup <- function(i) {
        print(i/ncells)
        tfile = paste(tfile0, 'analyseGroup', switch, model, Temp, i, ntest, '.Rd', sep = '-')
         
        if (file.exists(tfile)) {load(tfile); return(out)}
        ids = index[i:(min(i+ntest-1, length(index)))]
        out = sapply(ids,analyseCell, dat, 
                  gwt, Temp = Temp, ID =which.min(abs(gwt - Temp)))
        
        save(out, file = tfile)
        out
    }
    vout = lapply(seq(1, length(index), by = ntest), analyseGroup)
    voutAll = do.call(cbind, vout)
    out = dat[[1]][[1:4]]
    out[] = NaN
    fgwt = state = out[[1]]
    out[index] = t(voutAll)
    fgwt[] = gwt[out[[3]][]]
    fgwt[out[[3]][] == length(gwt)] = -1 

    above = out[[1]] > 1; below = !above; up = out[[2]]; down = !up
    state[above & up] = 3
    state[above & down] = 4
    state[below & up] = 2
    state[below & down] = 1
    
    outs = addLayer(out[[1]], state, fgwt)
    outs = writeRaster(outs, file = tfile)
    return(outs)
}

forTemp <- function(Temp, ...)
    out = mapply(analyseModel, dats, models, Temp, ...)

Temps = c(1, 1.5, 2)
outs = lapply(Temps, forTemp)
outsF = lapply(Temps, forTemp, switch = TRUE)


##########
## plot ##
##########
colsI = c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#f7f7f7','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')
limsI = c(0.1,  0.2, 0.4, 0.6, 0.8, 1, 1.25, 1.6, 2.5, 5, 10)
colsS = c("#8c510a", "#c7eae5", "#01665e", "#f6e8c3")
nlimT = 8
colsT = list(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7'),
             c('#d1e5f0','#92c5de','#4393c3','#2166ac','#053061'))

out = outs[[3]]
Temp = 2


#a = -log(0.5)/Temp
#limsT = -log(1-seq(0.1, 0.9, 0.1))/a
#limsT[limsT < 2] = round(limsT[limsT <2], 1)
#limsT[limsT > 2 & limsT < 3] = 0.5*round(2*limsT[limsT >2 & limsT < 3], 0)

#limsT[limsT > 3] = round(limsT[limsT >3], 0)
limsT = c(-2, -1, -0.5, -0.2, -0.1, 0, 0.1, 0.2, 0.5, 1, 2)
maxT = min(apply(gwts[,grepl('6.0', names(gwts))], 2, max))

limsT = limsT[limsT <= (maxT-Temp)]
colsT = c(make_col_vector(colsT[[1]], limits = limsT[limsT <=0]),
          make_col_vector(colsT[[2]], limits = limsT[limsT >=0]))
plotForTemp <- function(Temp, out, outF) {
png(paste0("figs/degreeEquil-", Temp, ".png"), height = 20, width = 12, units = 'in', res = 300)
layout(rbind(1:4, 5:8, 9, 10:13, 14:17, 18, 18+matrix(1:20, ncol = 4), 39, 40:43, 44), 
       heights = c(1, 1, 0.3, 1, 1, 0.3, rep(1, 5), 0.3, 1, 0.3))
par(mar = rep(0, 4), oma = c(0, 2, 2, 0))

plotImpact <- function(res, name, mttext = "Impact", addModName = TRUE){
    plotStandardMap(res[[1]], colsI, limsI)
    if (addModName) mtext(name, side = 3, line = 0)
    if (name == models[1]) mtext(mttext, side = 2)
}

mapply(plotImpact, out, models, "Impact fire off")
mapply(plotImpact, outF, models, "Impact fire on", FALSE)
StandardLegend(colsI, limsI, out[[1]][[1]])

plotState <- function(res, name, mttext) {
    plotStandardMap(res[[2]], colsS, seq(1.5, 3.5))
    if (name == models[1]) mtext(mttext, side = 2)
}

mapply(plotState, out, models, "Fire off")
mapply(plotState, outF, models, "Fire on" )

plot.new()
legend('center', horiz = TRUE, pch = 19, col = colsS, 
      c('reducing', 'recovering', 'increasing', 'diminishing'), bty = 'n')


plotTemp <- function(res, name) {  
    
    gwt = res[[3]] 
    mask = gwt == -1
    gwt = gwt - Temp

    state = res[[2]]

    plotTempState <- function(i) {
        if (i > 0) {
            if (sum(state[] == i, na.rm = TRUE) == 0) {plot.new(); return()}
            gwt[state != i] = NaN
        }
        
        plotStandardMap(gwt, colsT, limsT)
        plotStandardMap(mask, cols = c("transparent", "#9999aa"), limits = c(0.5), add = TRUE)
          
        if (name == models[1]) mtext(c('All', 'reducing', 'recovering', 
                                       'increasing', 'diminishing')[i+1], side = 2)
    }
    lapply(0:4, plotTempState)
}

mapply( plotTemp, out, models)
StandardLegend(colsT, limsT, out[[1]][[1]], extend_min = TRUE)
dev.off()
}

mapply(plotForTemp, Temps, outs, outsF)