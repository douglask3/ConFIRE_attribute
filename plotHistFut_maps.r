###########
## setup ##
###########
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/")
library(ncdf4)
sourceAllLibs("libs/")
graphics.off()

############
## params ##
############
dir = '../ConFIRE_ISIMIP/outputs/sampled_posterior_ConFire_ISIMIP_solutions/attempt4-full/'
models = c("GFDL-ESM2M", "HADGEM2-ES", "MIROC5", "IPSL-CM5A-LR")
periods = c("historic", "RCP2.6_2090s", "RCP6.0_2090s")
variable = "burnt_area_mean"

limits =  c(0,1, 2, 5, 10, 20, 40)
dlimits = c(-10, -5, -2, -1, -0.1, 0.1,  1, 2, 5, 10)
sdlimits =  c(-95, -90, -75, -50, 50, 75, 90, 95)
cols = c('#ffffcc','#ffeda0','#fed976','#feb24c','#fd8d3c','#fc4e2a','#e31a1c','#bd0026','#800026')

dcols = rev(c('#a50026','#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4','#313695'))

sc = 12 * 100

obs_file = "../jules_benchmarking/data/GFED4s_burnt_area.nc"


##########
## Open ##
##########
egRaster = raster('../ConFIRE_ISIMIP/outputs/sampled_posterior_ConFire_ISIMIP_solutions/attempt4-full/HADGEM2-ES/historic/sample_no_0.nc')
obs = raster::resample(mean(brick(obs_file)), egRaster)
obsL = (obs)

OpenPlotMap <- function(model, period, cols, limits, dcols = NULL, dlimits = NULL, dat0 = NULL,
                    anomolise = NULL, pnew = TRUE) {
    file = paste(dir, model, period, "fullPost.nc", sep = '/')
    file = paste(dir, model, period, "model_summary-50.nc", sep = '/')
    if (!file.exists(file)) file = paste(dir, model, period, "model_summary-10.nc", sep = '/')
    dat = brick(file, varname = variable)[[c(16, 50, 84)]]
    
    if (!is.null(anomolise)) 
        if (is.null(dat0)) datP = obs else {        
            dat = logistic(logit(dat) - anomolise)
            datP = dat
    } else datP = dat
    plotStandardMap(datP * sc, cols = cols,limits = limits)
    if (model == models[1]) mtext(side = 3, period)
    if (is.null(dat0)) {
        mtext(side = 2, adj = 0.2, model)
        if (pnew) plot.new()
    } else {
        dat = (dat-dat0) * sc
        plotStandardMap(dat, cols = dcols, limits = dlimits)
    }
    return(dat)
}

##########
## plot ##
##########
plotFun <- function(fname, anomolise = FALSE, signify = TRUE, controlT = TRUE) {
    print(fname)
    print(variable)
    print("---")
    #######################
    ## open plot all map ##
    #######################
    pdf(paste0("figs/", fname, "anomIs_", anomolise, "_Ymaps.pdf"),
        height = 10, width = 7.2)#, units = 'in',res  = 300)
        layout(rbind(matrix(1:(6*length(models)), ncol = 3), (6*length(models)) + c(1, 2, 2)))
        par(mar = rep(0, 4), oma = c(0, 2, 2, 0))
        
        if (anomolise) datA = obs else datA = NULL
        dat0 = lapply(models, OpenPlotMap, periods[1], cols, limits, anomolise = datA)
        
        if (anomolise) {
            datA = lapply(dat0, function(i) logit(i) - logit(obs))
            dat0 = list(obs)
        } else datA = list(NULL)
        
       
        datP = lapply(periods[2:3], function(p)
                    mapply(OpenPlotMap, models, dat0 = dat0, anomolise = datA,
                    MoreArgs = list(p, cols, limits, dcols, dlimits)))
        StandardLegend( cols,  limits, dat0[[1]])
        StandardLegend(dcols, dlimits, datP[[1]][[1]],
                            extend_min = TRUE, extend_max = TRUE)       
            
    dev.off.gitWatermark()
    pvsp = list(NULL)

    #############################################
    ## Signicance of difference from historic ##
    ############################################
    if (signify) {        
        signifM <- function(mid) {
            if (controlT) FUN <- function(x) x else FUN = logit
            dx = 0.01
            sce = 0.68
            xs = seq(-23, 23, dx)
            h = dat0[[mid]][[2]]
            #h = FUN(dat0[[mid]])[[2]]
            mask = !(h==h[1] & is.na(obsL))
            if (controlT) {
                error = sd((dat0[[mid]][[3]]-dat0[[mid]][[1]])[],
                          na.rm = TRUE )
            } else {
                error = sd((FUN(h) -obsL)[], na.rm = TRUE)*sce
            }
            h[!mask] = NaN          
            out = h
            #h = addLayer(h - error*1.5, h + error*1.5)
            signifP <- function(pid) {
                #### This bit changes when full model post is done 
                temp_file = paste("temp/pvs3", mid,  pid, sce,anomolise,fname,'.nc', sep = '-') 
                print(temp_file) 
                if (file.exists(temp_file)) return(raster(temp_file))
                #ft = FUN((datP[[pid]][[mid]][[2]]+dat0[[mid]][[2]]*sc)/sc)
                
                #ft = addLayer(ft - error*1.5, ft + error*1.5)
                #FUN4cell <- function(x, y) {
                #    p = dnorm(xs, x, error)
                #    q = dnorm(xs, y, error)
                #    sum(p*log(p/q))
                #}
                
                #pv = 1-exp(-mapply(FUN4cell, h[mask], ft[mask]))
                #out[mask] = pv
                #test = datP[[pid]][[mid]][[2]] < 0
                #out[test] = -out[test] 
               
                out = cal_pd_product(xs, h, datP[[pid]][[mid]][[2]],
                                     sc, error, mask, FUN)        
                #browser()
                out = writeRaster(out, file = temp_file, overwrite = TRUE)   
                        
            }
            lapply(1:2, signifP)    
        }
        pvs = lapply(1:length(dat0), signifM)
        #sigsP = lapply(pvs, lapply, function(i) i>  0.9 )
        #sigsN = lapply(pvs, lapply, function(i) i<(-0.9))
        #browser()

        out = pvs[[1]][[1]]
        out[] = NaN

        ## Under fire
        Sfire = max(layer.apply(pvs, layer.apply, function(i) abs(i))) > 0.9
        out[Sfire] = 1

        # under climate        
        RCPs = lapply(1:2, function(i) mean(layer.apply(pvs, function(p) p[[i]])))
        Smods = (abs(RCPs[[1]])>0.9) | (abs(RCPs[[2]]) > 0.9)
        out[Smods] = 2

        mods = layer.apply(pvs, function(i) abs(mean(layer.apply(i, function(x) x))))
        Srcp = max(mods) > 0.9
        out[Srcp] = 3
               
        All = mean(layer.apply(pvs, layer.apply, function(i) i))
        AllP = All > 0.9
        AllN = All < -0.9
        out[AllP] = 4
        out[AllN] = 5

        
        #browser()
        ## All agree
        #for (md in 1:length(sigsP)) for (rcp in 1:length(sigsP[[1]]))
        png("figs/Siggys.png", height = 5, width = 7.2, res = 300, units = 'in')  
        colsSigs = c('#fc8d62', '#8da0cb', '#66c2a5', '#a50026', '#313695')
        plotStandardMap(out, readyCut = TRUE, cols = colsSigs,
                        limits = seq(1.5, 4.5))
        legend('bottomleft', col = colsSigs, pch = 19, pt.cex = 2, legend = c('Sig. fire +', '\t climate', '\t RCP', '\t All - Increase', '\t All - Decrease'), bty = 'n')
        dev.off()

        out[] = NaN
        out[RCPs[[2]] >  0.9 & RCPs[[1]] <  0.9] = 1
        out[RCPs[[2]] <  0.9 & RCPs[[1]] >  0.9] = 2
        out[RCPs[[2]] < -0.9 & RCPs[[1]] > -0.9] = 3
        out[RCPs[[2]] > -0.9 & RCPs[[1]] < -0.9] = 4
     
        png("figs/Mitigation.png", height = 5, width = 7.2, res = 300, units = 'in')  
        colsSigs = c('#ca0020', '#92c5de', '#f4a582', '#0571b0')
        plotStandardMap(out, readyCut = TRUE, cols = colsSigs,
                        limits = seq(1.5, 4.5))
        legend('bottomleft', col = colsSigs, pch = 19, pt.cex = 2, legend = c('Sig. fire +', '\t climate', '\t RCP', '\t All - Increase', '\t All - Decrease'), bty = 'n')
        dev.off()

        browser()         
#browser() 
        fname0  = fname
        fname = paste0(fname, "anomIs_", anomolise, "signifChange", signify, "_Ymaps")
        pdf(paste0("figs/", fname, "pdf"), height = 5, width = 7.2)#, units = 'in',res  = 300)
         

            layout(cbind(1:5, c(6:9, 14), c(10:13,  14)), height = c(1,1, 1, 1, 0.5))
            par(mar = rep(0, 4), oma = c(0, 2, 2, 0))
            lapply(models, OpenPlotMap, periods[1], cols, limits,
                   anomolise = datA, pnew = FALSE)
            
            StanrdardLegend.new( cols,  limits, dat0[[1]])
            pvsp = lapply(1:length(pvs[[1]]), function(i) lapply(pvs, function(j) 100*j[[i]]))
           
            lapply(pvsp, function(i)
                        lapply(i, plotStandardMap, cols = dcols, limits = sdlimits))           
            mtext(outer = TRUE, "RCP2.6", adj = 0.5)
            mtext(outer = TRUE, "RCP8.0", adj = 0.85)   
            StanrdardLegend.new(dcols, sdlimits, pvsp[[1]][[1]], 
                                extend_min = TRUE, extend_max = TRUE)
        dev.off()
    }  
    
    triangulaise <- function(p, addLegend, lab, dlimits, pval = NULL) {
        perctile <- function(i) {            
            if (nlayers(p[[1]])>1) q = layer.apply(p, function(q) q[[i]]) else q = p
            P = mean(layer.apply(q, function(i) {i[i<0] = 0; i}))
            N = abs(mean(layer.apply(q, function(i) {i[i>0] = 0; i})))
            fout = paste0("outputs/", fnameP,'_',lab, '.nc')
            writeRaster.gitInfo(addLayer(P, N), file = fout, overwrite = TRUE)
            
            PNcols = make_col_vector(dcols, ncols = length( dlimits)+1)
            
            PNlims = dlimits[dlimits >0]
            P = cut_results(P, PNlims)
            N = cut_results(N, PNlims)
            PN = P + (N - 1) * (length(PNlims) + 1)
            
            k = 0
            PNcols = c()
            cols = make_col_vector(c(Ncols, Pcols), ncols = 101)
            index = seq(0, 1, length.out = length(PNlims)+1)
            dleg = (index[2]- index[1])/2
            legXY = c()
            for (i in index) for (j in index) {
                if (i == 0 && j == 0) rt = 50 else rt = round(100*j/(j+i))
                
                col = col0 = cols[rt+1]
                mag = round(100*(sqrt(i^2 + j^2)^0.5))
                mag[mag > 100] = 100
                col = make_col_vector(c("white", col), ncols = 101)[mag+1]
                #if (is.na(col)) browser()
                PNcols = c(PNcols, col)
                legXY = rbind(legXY, c(i, j))
            }            
            if (is.null(pval) || !controlT) e = NULL
            else e = 4 - sum(layer.apply(pval, function(i) abs(i)>95))
            
            plotStandardMap(PN, e = e, limits_error = 0.5:3.5, cols = PNcols,
                            limits = 0.5+1:((length(PNlims)+1)^2-1), ePatternRes = 50, 
                            ePatternThick = 0.67,
                            readyCut = TRUE, speedy = FALSE)
            mtext(lab, side = 2, line = -8)
            if (addLegend) {
                
                par(mar = c(3, 10, 1, 33.5))
                ext = range(index) + c(-dleg, dleg*1.1)
                plot(ext, ext, axes = FALSE, xlab = '', ylab ='', xaxs = 'i', yaxs = 'i')
                addSquare <- function(xy) {
                    if (as.numeric(xy[2]) < (1+1.5*dleg-as.numeric(xy[1]))) 
                        polygon(as.numeric(xy[1]) + dleg * c(-1, 1, 1, -1),
                                as.numeric(xy[2]) + dleg * c(-1, -1, 1, 1), col = xy[3])
                }
                apply(cbind(legXY, PNcols), 1, addSquare)
                #for (cex in 10*seq(10, 0, by = -0.2))
                #    points(legXY[,1], legXY[,2], pch = 19, cex = cex, col = PNcols)
                
                #polygon(c(0, 1, 1, 0), c(1, 1, 0, 1), col = "white", border = "white")
                mtext("Decrease (%)", side = 1, line = 2)
                mtext("Increase (%)", side = 2, line = 2)
                
                at1 = index - dleg#c(index[1] - diff(index[1])/2, index + diff(index)/2)
                at2 = rep(-0.12-dleg, length(at1))
                labels = c(0, PNlims)
                
                axisFUN <- function(id) {
                    #text(x = index[1], y = index[1], '0', xpd = NA, cex = 0.9)
                    text(x = at1, y = at2, labels, srt = 45, xpd = NA, cex = 0.9)
                    text(y = at1, x = at2, labels, srt = 45, xpd = NA, cex = 0.9)
                    lapply(1:2, axis, at = (index - diff(index)/2),
                                    labels = rep('', length(index)))
                    #id = seq(id, length(labels), by = 2)
                    #lapply(1:2, axis, at = (index - c(0, diff(index)/2))[id],
                    #                labels = labels[id], srt = 45)
                }
                lapply(1:2, axisFUN)
                lapply(1:2, axis, at = c(-9e9, 9e9), labels = c('', ''))

            
            }
    
        }
        #layer.apply(1:3, perctile)
        perctile(2)
    }
    if (controlT) {
        datP = lapply(datP, function(rs) mapply(function(r, r0)
                            r/max.raster(r0[[1]], na.rm = T), rs, dat0))
        dat0 = lapply(dat0, function(r) r/max.raster(r[[1]], na.rm = T))
    }
    fnameP = paste0(fname, "anomIs_", anomolise, "_CXmaps")
    pdf(paste0("figs/", fnameP, ".pdf"),
        height = 10, width = 7.2)#, units = 'in',res  = 300)
        layout(rbind(1, 2, 3, 4, c(5, 4), c(5, 5), 6), heights = c(1, 0.3, 1, 0.58, 0.42, 0.44, 0.56), widths = c(0.25, 0.75))
        par(mar = rep(0, 4), oma = c(0.1, 0, 0, 0))
        if (length(dat0) > 1)
            dat0 = mean(layer.apply(dat0, function(i) i[[2]]))
        else dat0 = dat0[[1]]
        writeRaster.gitInfo(dat0 * sc, file = paste0('outputs/', fnameP, '-GFED4s', '.nc'),
                            overwrite = TRUE)
        plotStandardMap(dat0 * sc, cols = cols, limits = limits, speedy = FALSE)
        mtext(side = 2, "Historic", line = -8)
        if (controlT) {
            extend_max = FALSE
            maxLab = 100
        } else {
            extend_max = TRUE
            maxLab = NULL
        }
        
        StanrdardLegend.new(cols,  limits, dat0[[1]], extend_max = extend_max, maxLab = maxLab, 
                  units = '%')
        
        mapply(triangulaise, datP, pval = pvsp, c(F, T), c("RCP2.6", "RCP6.0"),
              MoreArgs = list(dlimits = dlimits))
    dev.off()
   
    if (signify) {
        fnameP = paste0(fname, "anomIs_", anomolise, "_CXmaps")
        pdf(paste0("figs/", fnameP, "anomIs_", anomolise, "sign",signify, "_CXmaps.pdf"),
            height = 10, width = 7.2)#, units = 'in',res  = 300)#height = 7, width = 7.2)#, units = 'in',res  = 300)
            #par(mar = rep(0, 4), mfrow = c(3, 1))
            layout(rbind(1, 2, 3, 4, c(5, 4), c(5, 5), 6), heights = c(1, 0.3, 1, 0.58, 0.42, 0.38, 0.62), widths = c(0.25, 0.75))
            par(mar = rep(0, 4), oma = c(0.1, 0, 0, 0))
            #browser()
            plot.new()
            plot.new()
            mapply(triangulaise, pvsp, c(F, T), c("RCP2.6", "RCP6.0"),
                   MoreArgs = list(dlimits =sdlimits))
        dev.off()
        
    }  
}           
Pcols = "#330000"
Ncols = "#000033"

plotFun("BA", controlT = FALSE)

limits =  c(0, 0.1, 0.5, 1, 5, 10, 50)
dlimits = c(-10, -5, -1, -0.5, -0.1, -0.05, -0.01, 0.01, 0.05, 0.1, 0.5, 1, 5, 10)

plotFun("BA", TRUE, FALSE)

plotCOntrols <- function(control, 
                         slimits = c(0, 1, 2, 5, 10, 20, 50, 100, 150)/10,
                        sdlimits = c(-50, -20, -10, -5, -2, -1, 1, 2, 5, 10, 20, 50)/10) {
    variable <<- paste0("standard_",control)
    sc <<-  100# * 10
    limits <<-slimits
    dlimits <<- sdlimits
    plotFun(variable)

    sc <<-  100 * 12
    variable <<- paste0("potential_",control)
    limits <<- c(0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.15)
    dlimits <<-  c(-0.2, -0.1, -0.01, -0.001, 0.001, 0.01, 0.1, 0.2)/10
    #plotFun(variable)

    variable <<- paste0("sensitivity_",control)
    
    sc <<- 100 * 4 * 12
    #limits <<- c(0, 1, 2, 5, 10, 20, 50, 100, 150)/10
    #limits <<- c(-100, -50, -20, -10, -5, -2, -1, 1, 2, 5, 10, 20, 50, 100)/10
    if (control == "moisture") variable <<- "unknown"
    #plotFun(variable)    
}

cols = c('#ffffe5','#f7fcb9','#d9f0a3','#addd8e','#78c679','#41ab5d','#238443','#006837','#004529')
dcols = c('#40004b','#762a83','#9970ab','#c2a5cf','#e7d4e8','#f7f7f7','#d9f0d3','#a6dba0','#5aae61','#1b7837','#00441b')

Ncols = "#200026"
Pcols = "#002226"

plotCOntrols("fuel", slimits = c(5, 10, 20, 50, 80, 90, 95),
            sdlimits = c(-20, -10, -5, -2, -1, 1, 2, 5, 10, 20))

cols = c('#ffffe5','#fff7bc','#fee391','#fec44f','#fe9929','#ec7014','#cc4c02','#993404','#662506', '#331303')

dcols = rev(c('#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30'))

Pcols = "#131600"
Ncols = "#001620"

plotCOntrols("moisture", slimits = c(5, 10, 20, 50, 80, 90, 95),
            sdlimits = c(-20, -10, -5, -2, -1, 1, 2, 5, 10, 20))
