source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs('../rasterextrafuns/rasterPlotFunctions/R/')
source("libs/return_multiple_from_functions.r")

library(plotrix)
library(mapdata)
library(mapplots)

StandardLegend <- function(cols, limits, dat, rightx = 0.95, extend_max = TRUE, ...) 
        add_raster_legend2(cols, limits, dat = dat, add = FALSE,
                           transpose = FALSE, srt = 0, oneSideLabels= TRUE,
                           plot_loc = c(0.01, rightx, 0.3, 0.78),
                           ylabposScling = 1, extend_max = extend_max, ...)

plotStandardMap <- function(r, cols, limits, e = NULL, add_legend = FALSE,
                            limits_error = c(0.05, 0.1),
                            title2 = '', title3 = '', ...) {
    if (nlayers(r) > 1 && is.null(e)) {
        e = sd.raster(r)
        r = mean(r)
    } 
    r[r>9E9] = NaN
    if (!is.null(e)) e[is.na(r)] = NaN
    plot_raster_from_raster(r, e = e,
                            cols = cols, limits = limits, add_legend = FALSE,
                            quick = TRUE, ePatternRes = 5, ePatternThick = 0.67,
                            limits_error = limits_error, ...)
                            
    polygon(c(-62.5, -35, -35, -62.5), c(-56, -56, -50, -50), border = NA, col = "white")
    mtext(title3, adj = 0.1, line = 0.0)
    mtext(title2, side = 2, line = -1.5)
    if (add_legend) {
        add_raster_legend2(cols, limits, dat = trend[[2]],
                           transpose = FALSE, srt = 0, oneSideLabels= FALSE,
                           plot_loc = c(0.35, 0.99, 0.09, 0.12),  ylabposScling=0.8, ...)
    }
}