dat = read.csv("data/Tas_vs_year-ISIMIP.csv", stringsAsFactors=FALSE)[,1:5]


y_range = range(dat[,-1])
x_range = range(dat[,1])

png("figs/globaTemp-isimip2b.png", width = 7, height = 5.5, units = 'in', res = 300)
plot(x_range, y_range, type = 'n', axes = F, xlab = '', ylab = '')

grid()

cols = c('black', '#1b9e77','#d95f02','#7570b3')
plotDat <- function(i, ...) {
    lines(dat[,1], dat[,i+1], col = cols[i], lwd = 2, ...)
}

lapply(1:4, plotDat)
lapply(1:4, plotDat, lty = 2)

modnames = colnames(dat)[-1]
modnames = substr(modnames, 1, nchar(modnames) - 4)
legend('topleft', lty = 1, col = cols, modnames, lwd = 2)

source("../rasterextrafuns/rasterPlotFunctions/R/mtext.units.r")
axis(1)
axis(1, at = c(1850, 2100))
axis(2)
axis(2, at = c(-9, 9))
mtext.units(side = 2, 'GWL(~DEG~C)', line = 2)
dev.off()
