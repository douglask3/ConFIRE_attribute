################################################################################
## cfg                                                                        ##
################################################################################
## Libraries etc
#source('cfg.r')
library(rstash)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
 source("libs/citation.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")

## paths and parameters
dir   = 'data/climate/cru_ts4.05/'
varns = c(temp   = 'tmp',
          precip = 'pre',
          cloud  = 'cld')

clim_layers = 949:1440
################################################################################
## load data                                                                  ##
################################################################################
files = list.files(dir, full.names = TRUE)

clim_layers= (min(clim_layers-12):max(clim_layers))
loadDat <- function(varn) {
    files = files[grepl(varn, files)]    
    return(layer.apply(files, stack)[[clim_layers]])
}

dat = lapply(varns, loadDat)

nyears = floor(length(clim_layers)/12)

################################################################################
## Function for stash                                                         ##
################################################################################
## convert data to xyz format ready for stash
r2xyz <- function(dat, index) {
    dat = dat[[index]]
    xy  = xyFromCell(dat[[1]], 1:length(dat[[1]]))
    dat = cbind(xy, values(dat))
    return(data.frame(dat))
}

## run stash for the year
year.stash <- function(y, spinup = FALSE) {
    cat('Year ', y, 'of ', nyears, '\n')

    # select and xyz required timesetps
    m = (1 + (y - 1) * 12):(y*12)
    dat0 = dat[[1]][[m]]
    dat = lapply(dat, r2xyz, m)

    # set up cell data
    gchar = dat[[1]][1:2]
    gchar$elev = 0
    gchar$fcap = 140
    # initial soil water conditions [set -9999.0 for spinup]
    if (exists('swc0')) gchar$swc0 = swc0
        else gchar$swc0 = -9999.0

    stash = grid.stash(dat[[1]], dat[[2]], 1 - dat[[3]]/100, gchar)
    alpha = stash$alpha.index
    swc0 <<- stash$swc.init[,3]

    if (!spinup) {
        dat0[] = as.matrix(alpha[,c(-1,-2)])
        return(dat0)
    }
    invisable()
}

################################################################################
## run and output                                                             ##
################################################################################

lapply(rep(1, 40), year.stash)
alpha = layer.apply(1:nyears, year.stash)
nms = names(dat[[1]])
names(alpha) = nms
date = as.Date(gsub('.', '-', substr(nms, 2, 11), fixed = TRUE))
alpha = setZ(alpha, date, 'Date')


#alphaMax = seaCy12(alpha, function(...) max(...))

comment = list('made using rstash' = citation.text('rstash'))
writeRaster.gitInfo.time <- function(...) 
    writeRaster.gitInfo(..., zunit = 'month', zname = 'time')
#alpha2000_2010 = alpha[[-c(1:6)]]
#alpha2000_2010 = alpha2000_2010[[1:132]]
writeRaster.gitInfo.time(alpha[[13:nlayers(alpha)]], "outputs/alpha1980-2020.nc",
                    comment = comment, overwrite = TRUE)

browser()
alpha = alpha[[-c(1:12)]]



writeRaster.gitInfo.time(alpha, drive_fname['alpha'],
                    comment = comment, overwrite = TRUE)

writeRaster.gitInfo.time(alphaMax, drive_fname['alphaMax'],
                    comment = comment, overwrite = TRUE)

writeRaster.gitInfo.time(alpha, drive_fname['alpha'],
                    comment = comment, overwrite = TRUE)

