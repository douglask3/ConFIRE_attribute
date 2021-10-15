StanrdardLegend.new <- function(cols, limits, dat, ...) 
    add_raster_legend2(cols, limits, dat = dat,
                           transpose = FALSE, srt = 0, oneSideLabels= TRUE,
                           plot_loc = c(0.1, 0.9, 0.73, 0.8), ylabposScling=0.8,
                           add = FALSE, ...)

StandardLegend <-  function(cols, limits, dat, rightx = 0.95, extend_max = TRUE, 
                            oneSideLabels = TRUE, add = FALSE, 
                            ylabposScling = 1,
                            ytopScaling = -1,  ...) {
    if (add)        
        plot_loc = c(0.41, rightx, 0.1, 0.13)
    else 
        plot_loc = c(0.01, rightx, 0.3, 0.56)
    add_raster_legend2(cols, limits, dat = dat, add = add,
                       transpose = FALSE, srt = 0, oneSideLabels= oneSideLabels,
                       plot_loc = plot_loc, ytopScaling = -1, 
                       ylabposScling = ylabposScling, extend_max = extend_max, ...)
}
