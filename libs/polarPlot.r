source("libs/make_transparent.r")

polarPlot.setup <- function(x, y, xlim = NULL, ylim = xlim, type = 'l',...) {
    if (type == 'l') {
        x = c(x, x[1])
        y = c(y, y[1])
    }
    x = 2 * pi * x/12
    yn = y * cos(x)
    xn = y * sin(x)
    x = xn; y = yn
    if (is.null(xlim)) xlim = ylim = max(abs(x), abs(y)) * c(-1, 1)
    
    plot(x, y, pch = 19, cex = 2, xlab = '', ylab = '', axes = FALSE, xlim = xlim, ylim = ylim, type = type,...)
}

polarPlot.polygon <- function(x, y, col = "black", alpha = 0.67, border = TRUE, ...) {
    x0 = x; y0 = y
    
    x = 2 * pi * x/12
    yn = apply(y, 1, function(i) i * cos(x))
    xn = apply(y, 1, function(i) i * sin(x)) 
    
    x = rbind(xn, xn[1, ]); y = rbind(yn, yn[1, ])
    
    polygon(x, y, border = NA, col = make.transparent(col, alpha), ...)
    if (border) {
        lines(x[,1], y[,1], col = col)
        lines(x[,2], y[,2], col = col)
    }
}

polarPlot.addGuides <- function(xlim = c(-1, 1), ylim = xlim, axisMonth = 0) {
    lines(c(0, 0), ylim/2, lwd = 2)
    lines(ylim/2, c(0, 0), lwd = 2)
    
    
    at = seq(0, signif(xlim[2], 1), length.ou = 6)
    at = signif(at, 2)
    at = at[at <= xlim[2]]
    
    mnths = 2 * pi *((0.5:11.5)/12)
    xr = xlim[2] * sin(mnths) * 1.07
    yr = ylim[2] * cos(mnths) * 1.07
    text(x = xr, y = yr, c('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'))
    
    addRadin <- function(mnths, ...) {
        xr = xlim[2] * sin(mnths) 
        yr = ylim[2] * cos(mnths) 
        mapply(function(i,j) lines(c(0, i), c(0, j), col = make.transparent("white", 0.67), ...), xr, yr)
        mapply(function(i,j) lines(c(0, i), c(0, j), col = make.transparent("black", 0.67), lty = 2, ...), xr, yr)
    }
    addRadin(2 * pi *((0:12)/12))
    for (i in 1:4) addRadin(2 * pi *(c(2, 5, 8, 11)/12), lwd = 2)
    
    axisMonth = 2 * pi * (axisMonth + 0.5)/12
    addCirclegrid <- function(r) {
        xr = r * sin(seq(0, 2*pi, 0.01) + axisMonth)
        yr = r * cos(seq(0, 2*pi, 0.01) + axisMonth)   
        lines(xr, yr, col =  make.transparent("white", 0.33))
        lines(xr, yr, lty = 2,  col =  make.transparent("black", 0.33))
        if (r == 0) cex = 2 else cex = 4
        points(xr[1], yr[1], pch = 19, cex = cex, col = "white")
        text(y = yr[1], x = xr[1], r)                
    }
    lapply(at, addCirclegrid)
}