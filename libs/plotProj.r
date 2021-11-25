source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs('../rasterextrafuns/rasterPlotFunctions/R/')
source("libs/return_multiple_from_functions.r")

library(plotrix)
library(mapdata)
library(mapplots)

library(rgdal)


plotProj <- function(r, cols, limits, e = NULL, add_legend = FALSE,
                            limits_error = c(0.2, 0.200000001),
                            title2 = '', title3 = '', xlim = c(-180, 180), ylim = c(-60, 90),
                            ePatternRes = 67, ePatternThick = 0.5, greyBck = FALSE, 
                            bckCol = 'black', txtCol = 'white', extend_min = FALSE,
                            ...,  speedy = T, id = '', units = '', oneSideLabels = FALSE,
                            noCut = FALSE,
                            start = c(10, -120, 0), end  = c(10, 239, 0),
                            nFrames = round(360*2/3), rotBase = 5, projection = 'orthographic',
                            prj_parameters = NULL, frameSpeed = 9,
                            spt.cex_data = 0.2, spt.cex_ocean = 0.3, 
                            spt.cex_ice = spt.cex_ocean, spt.cex_bare = spt.cex_ocean) {
    
    print(id)
    r0 = r
    selectR <- function(r, mn, bs) {
         
        if (is.list(r)) return(selectR(r[[sample(1:length(r), 1)]], mn, bs))
        if (nlayers(r) == 1) {
            rn = r[[1]]
            
        } else if (nlayers(r) == 3) {
        
            if (min.raster(r, na.rm = TRUE) < 0) {
                print("anaomolie")
                ntrans = TRUE
                r = (r+100)/2
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
            lyr = (mn %/% bs) %% 12 + 1
            rn = r[[lyr]]
           
        } else if(nlayers(r) > 1) {
            lyr = (mn %/%  bs)
            if (bs>1) lyr = lyr + 1
            rn = r[[lyr]]
            print(lyr)
        }
        
        if (!noCut) {
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
        } else {
            rn[rn == 0] = NaN
        }
        
        return(rn)
    }
    
    bareCols = c("#00441b", "#7fbc41", "yellow", "brown")
    dir = paste0('figs/', id)
    makeDir(dir)

    if (greyBck) bareCols = c('grey', 'white')
    
    plotFUN <- function(rot, ...) {
        #browser()
        nrots = paste(c(rep('0', 6-nchar(rot[1])), rot[1]), collapse = '')
        figname = paste0(dir, '/map-', nrots, '.png')
        if (file.exists(figname)) return()
        png(figname, height = 6, width = 6, units = 'in', res = 150)
        
        par(mar = rep(5, 4))
        FUN <- function(x, colsi = cols, limitsi = limits, add = TRUE, spt.cex=1){
            #browser() 
            plot_raster_from_raster(x,coast.lwd = 0, e = NULL,
                                    cols = colsi, limits = limitsi, add_legend = FALSE,
                                    quick = TRUE, add = add,spt.cex= spt.cex, 
                                    orientation = rot[-1], ...)
        }
        if (rot[1] == 1 || rot[1] %% rotBase == 1 || rotBase == 1) rn = selectR(r, rot[1], rotBase)
        #browser()
        rn <<- rn
    
       
        FUN(rn, add = FALSE)
        
        if (!is.null(bckCol)) {
            polygon(9E9 * c(-1, 1, 1, -1), 9E9 * c(-1, -1, 1, 1), xpd = NA, col =  bckCol)
            FUN(rn)
        }
        FUN(bare, bareCols, 0:100)
        
        
        FUN(rn, spt.cex=spt.cex_data * 2)
        FUN(ice, c("transparent", "#EFEFFF"), c(0.5))
        FUN(sea, c("transparent", "blue"), limits = c(0.5))
        FUN(oceanDepth, c("black", "blue"), c(-9000, -5000, -3000, -2000, -1000), spt.cex=spt.cex_ocean)
        
        FUN(ice, c("transparent", "#EFEFFF"), c(0.5), spt.cex=spt.cex_ice)
        FUN(icesheet, c("white", "#FFAA88"), c(0, 500, 1000, 1500, 2000, 2500, 3000, 3500, 4000), spt.cex=spt.cex_ice/3)
        #browser()
        FUN(bare, bareCols, 0:100, spt.cex=spt.cex_bare)
        FUN(rn, spt.cex=spt.cex_data)
        mnthNames = c('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 
                      'Aug', 'Sep', 'Oct', 'Nov', 'Dec')
        mn = (rot[1] %/% rotBase) %% 12+1
        
        if (nlayers(r) == 12) 
            mtext(side = 1, line = -1.4, adj = 0.09, mnthNames[mn], cex = 2, col = txtCol)
        if (nlayers(r) > 12) {
            txt = tail(strsplit(strsplit(names(rn),'.', fixed = TRUE)[[1]][1], 'X')[[1]],1)
            mtext(side = 1, line = -1.4, adj = 0.09, txt, cex = 2, col = txtCol)
        }
        #browser()
        dev.off()
        #FUN()
    }
    #browser()
    rot = cbind(1:nFrames, mapply(seq, start, end, length.out = nFrames))
    dir = paste0(c(dir, '/', projection, '-', prj_parameters, '-', 
                   start, '-', end, '-', nFrames, '-', rotBase, '/'), collapse = '')
    makeDir(dir)
    
    if (F) {
    png(paste0(dir, '/legend.png'), height = 5.5, width = 7, res = 300, units = 'in')
        plot.new()
        polygon(9E9 * c(-1, -1, 1, 1), 9E9 * c(-1, 1, 1, -1), col = bckCol, xpd = NA)
        StandardLegend(dat = r[[1]], cols = cols, limits = limits, add = TRUE, 
                       ylabposScling = 0.3, fgCol = txtCol, bckCol = bckCol, 
                       oneSideLabels = oneSideLabels, units = units,
                       extend_min = extend_min)
    dev.off()  
    }
    #rot = seq(start, end, by = turnSize)
    figureName = paste0(dir, '/gmap-', frameSpeed, '.gif')
    if (file.exists(figureName)) return()
    
    apply(rot, 1, plotFUN, projection = projection, prj_parameters = prj_parameters)
    commd = paste0("convert -delay 1 ", dir, '/map*.png ', figureName)
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

