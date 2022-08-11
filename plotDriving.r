###########
## setup ##
###########
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
library(ncdf4)
sourceAllLibs("libs/")
graphics.off()

dir = "/data/users/dkelley/ConFIRE_ISIMIP/INFERNOonOff/"

regTest = "/scratch/cburton/ISIMIP_PAPER/KScubes/KSresults_PD//"#HadGEM2_NBP_TOE_FireOff.nc"

models = c("HadGEM2", "IPSL", "MIROC", "GFDL")

gwls = read.csv('data/Tas_vs_year-ISIMIP.csv')


regFiles = list.files(regTest, full.names = TRUE)
forModel <- function(mod) {
    files = regFiles[grepl(mod, regFiles) & grepl('30', regFiles)]
    
    openF <- function(test) {
        out = raster(files[test])
        out[out == 0] = NaN
        out[abs(out) >2089] = NaN
        out
    }
    noFire = openF( grepl('FireOff', files))
    Fire   = openF(!grepl('FireOff', files))
    browser()
    gwl = gwls[,grep(mod, colnames(gwls))[1]]
    gwl[length(gwl)] = NaN
    gwl[gwl < 0] = 0
    returnTemp <- function(r) {
        mask = !is.na(r)
        index = abs(r[mask])-1869
        
        index[index > length(gwl)] = length(gwl)
        index[index < 1] = 1
        r[mask] = gwl[index] * (r[mask] > 0)
        return(r)
    }
    gwlFon  = returnTemp(  Fire)
    gwlFoff = returnTemp(noFire)
    
    
    plotStandardMap(noFire , cols=colsYrs, limits = limitsYrs)
    mtext(side = 2, mod, line = 1, adj = 0.1)
    plotStandardMap(Fire   , cols=colsYrs, limits = limitsYrs)

    plotQuadMap <- function(r1, r2, cols1, cols2, limits) {
        browser()
    }
    #plotQuadMap(Fire, noFire, colsYrs, colsGwl)
    #browser()

    plotStandardMap(gwlFoff, cols=colsGwl, limits = limitsGwl)
    plotStandardMap(gwlFon , cols=colsGwl, limits = limitsGwl)
    
}

limitsYrs = c(2005, 2007, 2010, 2020, 2050, 2070)
limitsGwl = seq(0.2, 2.2, by = 0.4)
limitsYrs = c(rev(-limitsYrs), 2004, limitsYrs)
limitsGwl = c(rev(-limitsGwl), 0, limitsGwl)

colsYrs = c(rev(c('#492004', '#7f3b08','#b35806','#e08214','#fdb863','#fee0b6')),#'#f7f7f7',
            rev(c('#d8daeb','#b2abd2','#8073ac','#542788','#2d004b', '#1f002b')))
colsGwl = c(rev(c('#67001f','#b2182b','#d6604d','#f4a582','#fddbc7')),rev(c('#d1e5f0','#92c5de','#4393c3','#2166ac','#053061')))
png("figs/ToE_GWoE.png", height = 5, width = 8, res = 300, units = 'in')
    layout(rbind(t(matrix(1:16, nrow = 4)),c(17, 17, 18, 18)), heights = c(1, 1, 1, 1, 0.3)) 
    par(mar = rep(0, 4), oma = c(2, 2, 0, 0))
    lapply(models, forModel)
    StandardLegend(colsYrs, limitsYrs, c(0,1), extend_min = TRUE, oneSideLabels = FALSE)
    StandardLegend(colsGwl, limitsGwl, c(0,1), extend_min = TRUE, oneSideLabels = FALSE)
dev.off()
