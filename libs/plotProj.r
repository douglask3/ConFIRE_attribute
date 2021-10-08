source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs('../rasterextrafuns/rasterPlotFunctions/R/')
source("libs/return_multiple_from_functions.r")

library(plotrix)
library(mapdata)
library(mapplots)

library(rgdal)

#SA_ste <- readOGR(dsn = "data/South_America", layer = "South_America")
#rivers <- readOGR(dsn = "data/majorrivers_0_0", layer = "MajorRivers")

StandardLegend <- function(cols, limits, dat, rightx = 0.95, extend_max = TRUE, oneSideLabels = TRUE, add = FALSE, ...) {
    if (add)        
        plot_loc = c(0.41, rightx, 0.1, 0.13)
    else 
        plot_loc = c(0.01, rightx, 0.3, 0.56)
    add_raster_legend2(cols, limits, dat = dat, add = add,
                       transpose = FALSE, srt = 0, oneSideLabels= oneSideLabels,
                       plot_loc = plot_loc,
                       ylabposScling = 1, extend_max = extend_max, ...)
}

lineBox <- function(x, y, ...) 
    lines(c(x[1], x[2], x[2], x[1], x[1]), c(y[1], y[1], y[2], y[2], y[1]),...)

plotProj <- function(r, cols, limits, e = NULL, add_legend = FALSE,
                            limits_error = c(0.2, 0.200000001),
                            title2 = '', title3 = '', xlim = c(-180, 180), ylim = c(-60, 90),
                            ePatternRes = 67, ePatternThick = 0.5,
                            ...,  speedy = T, id = '') {
    
    print(id)
    r0 = r
    selectR <- function(mn, bs) {
        if (nlayers(r) == 1) {
            rn = r[[1]]
        } else if (nlayers(r) == 3) {
        
            if (min.raster(r, na.rm = TRUE) < 0) {
                print("anaomolie")
                ntrans = TRUE
                r = (r0+100)/2
            } else ntrans = FALSE
            r = logit(r/100)
            rn = r[[1]]
            selectVal <- function(x) {
                out = rnorm(1, x[2], abs(x[2]-x[1]))
                #if (x[3] < x[2]) out = -out
                if (is.na(out) ) browser()
                return(out)
            }
            test = r[[1]] < r[[2]]
            rn[test] = apply(r[test], 1, selectVal)
            
            rn0 = rn
            if (ntrans) rn =100*2*logistic(rn)-100
            else {
                rn = logistic(rn)*100
                rn[rn == min.raster(rn)] = NaN
            }
        } else if(nlayers(r) == 12) {
            browser()
        } else if(nlayers(r) > 1) {
            browser()
        }
        if (limits[1] < 0) {
            
            cutP = min(abs(limits))
            cutP =  max(cutP, quantile(abs(rn[]), 0.5, na.rm = TRUE))     
            rn[abs(rn) < cutP] = NaN
        } else {
            if(limits[1] == 0) cutP = limits[2]
            else cutP = limits[1]
            cutP = max(cutP, quantile(rn[], 0.5, na.rm = TRUE))             
            rn[rn < cutP] = NaN
        }
     
        
        return(rn)
    }
    
    bareCols = c("#00441b", "#7fbc41", "yellow", "brown")
    dir = paste0('figs/', id)
    makeDir(dir)
    plotFUN <- function(rot, ...) {
        #browser()
        nrots = paste(c(rep('0', 6-nchar(rot[1])), rot[1]), collapse = '')
        figname = paste0(dir, '/map-', nrots, '.png')
        if (file.exists(figname)) return()
        png(figname, height = 6, width = 6, units = 'in', res = 300)
        
        par(mar = rep(5, 4))
        FUN <- function(x, colsi = cols, limitsi = limits, add = TRUE, spt.cex=1) 
            plot_raster_from_raster(x,coast.lwd = 0, 
                                    cols = colsi, limits = limitsi, add_legend = FALSE,
                                    quick = TRUE, add = add,spt.cex= spt.cex, 
                                    orientation = rot[-1], ...)
        #browser()
        if (rot[1] %% rotBase == 0) rn = selectR(rot[1], rotBase)
        rn <<- rn
        FUN(rn, add = FALSE)
        FUN(bare, bareCols, 0:100)
        FUN(rn, spt.cex=0.4)
        FUN(ice, c("transparent", "#EFEFFF"), c(0.5))
        FUN(sea, c("transparent", "blue"), limits = c(0.5))
        FUN(oceanDepth, c("black", "blue"), c(-9000, -5000, -3000, -2000, -1000), spt.cex=0.3)
        
        FUN(ice, c("transparent", "#EFEFFF"), c(0.5), spt.cex=0.3)
        FUN(icesheet, c("white", "#FFAA88"), c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000), spt.cex=0.1)
        FUN(bare, bareCols, 0:100, spt.cex=0.3)
        FUN(rn, spt.cex=0.2)
        
        dev.off()
        #FUN()
    }
    start = c(10, -120, 0)
    end  = c(10, 239, 0)
    nFrames = 360
    projection = 'orthographic'
    prj_parameters = NULL
    rotBase = 5
    frameSpeed = 9
    rot = cbind(1:nFrames, mapply(seq, start, end, length.out = nFrames))
    dir = paste0(c(dir, '/', projection, '-', prj_parameters, '-', 
                   start, '-', end, '-', nFrames, '-', rotBase, '/'), collapse = '')
    makeDir(dir)
    
    #rot = seq(start, end, by = turnSize)
    figureName = paste0(dir, '/gmap-', frameSpeed, '.gif')
    if (file.exists(figureName)) return()
    
    apply(rot, 1, plotFUN, projection = projection, prj_parameters = prj_parameters)
    commd = paste0("convert -delay 1 ", dir, '/*.png ', figureName)
    system(commd)
   
}

addCoastlineAndIce2map <- function() {
    add_icemask()
    
    mask = raster('data/seamask.nc')
    mask = mask>1
    
    plot_raster_from_raster(mask+1, add = TRUE, 
                             cols = c("white", "transparent"),readyCut = TRUE,
                             limits =  NULL, quick = TRUE, interior = FALSE, 
                             coast.lwd = NULL, add_legend = FALSE)
    #
    #contour(mask, add = TRUE, drawlabels = FALSE, lwd = 0.5)  

    ployBox <- function(x, y)
        polygon(c(x[1], x[2], x[2], x[1]), c(y[1], y[1], y[2], y[2]), col = "white", border = "white")
        
    ployBox(c(-180, -90), c(-60, 0))
    ployBox(c(-180, -120), c(-60, 25))
    ployBox(c(-50, -19), c(10, 25))
    ployBox(c(-50, -13.5), c(27.5, 34))
    ployBox(c(115, 125), c(-8, -7))
    ployBox(c(104, 111), c(2.5, 8))
    ployBox(c(122, 128), c(2.5, 5)) 
}

add_icemask <- function() {
	icemask = raster('data/icemask.nc')
	plot_raster_from_raster(icemask, add = TRUE, cols = c('#FFFFFFFF', 'grey'), y_range = c(-60, 90),
						    limits = c(-0.5, 0.5), add_legend = FALSE, interior = FALSE, coast.lwd = 0.67)#, coast.lwd = NULL)
}

