graphics.off()
years = 1939:2089

temps = read.csv("data/Tas_vs_year-ISIMIP.csv")
temps = sapply(years, function(yr) temps[which.min(abs(temps[,1] -yr)),4])

yr1p5 = years[which.min(abs(temps - 1.5))]
wh1p5 = which(years == yr1p5)
ynorm = (years - min(years))/diff(range(years))

logitstic <- function(x) 1/(1-exp(-x))

newPlot <- function(case1, case2, case3, nme, letter = '') {
    #case2 = case2/case1[1]; case3 = case3/case1[1]; case1 = case1/case1[1]
    
    plot(range(years), c(0, 1), type = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
    axis(3)
    mtext(side = 3, letter, line = 0, adj = 0.01)
    temp_ticks = c(0, 0.5, 1, 1.5, 2, 2.5)
    tick_locs = sapply(temp_ticks, function(i) which.min(abs(temps-i)))
    axis(1, at = years[tick_locs], labels = temp_ticks)
    axis(2, at = case1[2], label = 1)
    axis(2, at = seq(0, 1, 0.2), labels = rep('', 6))
    lines(years, case1, col = 'blue')
    lines(years, case2, col = 'red')
    lines(years, case3, col = 'red')
    
    lines(years, case1, lty = 2, col = 'blue')

    lines(c(yr1p5, yr1p5), c(-10, case1[wh1p5]))
    lineSpan <- function(old, new) {
        newWh = old[wh1p5] - new
        testDir = diff(old[c(wh1p5-1, wh1p5)]) < 0
        possibles = c(which(newWh == 0 |
                            newWh[-1] < 0 & head(newWh, -1) > 0 | 
                            newWh[-1] > 0 & head(newWh, -1) < 0))

        testPoss <- function(i) (diff(new[c(i-1, i)]) < 0) == testDir
        possibles = possibles[sapply(possibles, testPoss)]
        
        newWh = possibles[which.min(abs(possibles - wh1p5))]
        lines(c(yr1p5, years[newWh]), rep(case1[wh1p5], 2), lty = 2)
        lines(rep(years[newWh], 2), c(case1[wh1p5], -0.2), lty = 2, xpd = TRUE)
    }

    lineSpan(case1, case2)
    lineSpan(case1, case3)
    mtext(side = 2, nme, line = 0.5)
}

png("tree_equl_explanation.png", height = 5, width = 7, res = 300, units = 'in')
par(mfrow = c(2, 2), mar = c(2, 1, 2, 1), oma = c(0.5, 1, 0.1, 0.2))
newPlot(1-ynorm+0.006666667 , 1-ynorm/1.5+0.006666667 , 1-ynorm*2+0.006666667 , 'Reducing', letter = 'a')
newPlot(ynorm-0.006666667 , ynorm/1.5-0.006666667 , ynorm*2-0.006666667 , 'Increasing', letter = 'b')

#yr1p5 = 2064
#wh1p5 = which(years == yr1p5)
#newPlot(((ynorm/0.5)-0.5)^2, ((ynorm/0.59)-0.8)^2, ((ynorm/0.35)-1)^2, 'Recovering')
newPlot(-3.4+2.7*(ynorm^(-0.1)+ ynorm - 0.2*(ynorm+0.5)^2)+0.0643397,
        -3.4+2.7*(ynorm^(-0.1)+ 0.86* ynorm - 0.2*(ynorm+0.5)^2)+0.0643397,
        -3.4+2.7*(ynorm^(-0.1)+ 1.2*ynorm - 0.2*(ynorm+0.5)^2)+0.0643397, 'Recovering', letter = 'c')

newPlot(4.4-2.7*(ynorm^(-0.1)+ ynorm - 0.2*(ynorm+0.5)^2) - 0.0643397,
        4.4-2.7*(ynorm^(-0.1)+ 0.86* ynorm - 0.2*(ynorm+0.5)^2) - 0.0643397,
        4.4-2.7*(ynorm^(-0.1)+ 1.2*ynorm - 0.2*(ynorm+0.5)^2) - 0.0643397, 'Diminishing', letter = 'd')

dev.off()
