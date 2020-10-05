library(ncdf4)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
source("libs/process_jules_file.r")

dirs = list(historic = '/hpc/data/d01/hadcam/jules_output/RUNC20C_u-bk886_isimip_0p5deg_origsoil_dailytrif/',
            RCP2.6 = '/hpc/data/d01/hadcam/jules_output/RUNFUT26_u-bk886_isimip_0p5deg_origsoil_dailytrif/',
            RCP6.0 = '/hpc/data/d01/hadcam/jules_output/RUNFUT60_u-bk886_isimip_0p5deg_origsoil_dailytrif/')

years = list(historic = 1995:2004,
             RCP2.6 = 2090:2099,
             RCP6.0 = 2090:2099)

fileIDs = c(cover = "ilamb", soilM = "gen_mon_layer", precip = "ilamb", humid = "ilamb", tas = "ilamb")

varnames =  c(cover = "frac", soilM = "smcl", precip = "precip", humid = "q1p5m_gb", tas = "t1p5m_gb")

models = c("GFDL-ESM2M", "HADGEM2-ES", "IPSL-CM5A-LR", "MIROC5")

temp_dir = '/data/users/dkelley/ConFIRE_ISIMIP_temp/-makeISIMIPins'
out_dir  = '/data/users/dkelley/ConFIRE_ISIMIP/inputs/'

coverTypes = list(trees = c(1:7), totalVeg = c(1:13), crop = c(10, 12), pas = c(11, 13))

makeDat <- function(id, dir, years_out) {
    years = c(years_out[1] - 1, years_out, tail(years_out, 1) + 1)
    forModel <- function(mod) {
    tfile0 = paste0(c(temp_dir, id, mod, range(years)), collapse = '-')
        dir = paste0(dir, '/', mod, '/')
        files = list.files(dir, full.names = TRUE)

        ## select years
        files = files[apply(sapply(years, function(i) grepl(i, files)), 1, any)]
    
        openVar <- function(fileID, vname) {
            tfile = paste(tfile0 , fileID, vname, '.Rd', sep = '-')
            
            if (file.exists(tfile)) {
                load(tfile)
            } else {
                print(tfile)
                files = files[grepl(fileID, files)]     
                dat = lapply(files, process.jules.file, NULL, vname)
                save(dat, file =  tfile)
            }
            return(dat)
        }
        dats = mapply(openVar, fileIDs, varnames)
        cover = dats[-1, 'cover']
        makeCover <- function(ty) {
            print(ty)
            group <- function(i) {
                cv = i[ty]
                out = cv[[1]]
                for (i in cv[-1]) out = out + i
                return(out)
            }
            coverTy = layer.apply(cover, group)
        }
        tfile = paste(tfile0, 'covers', '.Rd', sep = '-')
        if (file.exists(tfile)) load(tfile)
        else {
            covers = lapply(coverTypes, makeCover)
            save(covers, file = tfile)
        }
        tfile = paste(tfile0, 'soils', '.Rd', sep = '-')
        if (file.exists(tfile)) load(tfile)
        else {
            soilM_top    = layer.apply(dats[, 'soilM'], function(i) i[[1]])
            soilM_bottom = layer.apply(dats[-1, 'soilM'], function(i) i[[1]])

            st = 2:(nlayers(soilM_top)-11)
            ed = 13:nlayers(soilM_top)
            soil12 = mapply(function(s, e) sum(soilM_top[[s:e]]), st, ed)
            soil12 = layer.apply(soil12, function(i) i)
            soilM_top =  soilM_top[[-(1:12)]]
            soil12 = soilM_top/soil12
            save(soilM_top, soilM_bottom, soil12, file = tfile)
        }
        precip = layer.apply(dats[-1, 'precip'], function(i) i)
        tas = layer.apply(dats[-1, 'tas'], function(i) i)

        out_dirM = paste(out_dir,  mod, id, '', sep = '/')
        writeOut <- function(dat, name) {
            file = paste0(out_dirM,  name, '.nc')
            dat = dat[[-1]]
            nl = 12*floor(nlayers(dat)/12)
            dat = dat[[1:nl]]
           
            writeRaster.Standard(dat, file)
        }
        mapply(writeOut, covers, names(coverTypes))
        writeOut(soil12, 'soil12')
        writeOut(soilM_bottom, 'soilM_bottom')
        writeOut(soilM_top, 'soilM_top')
        writeOut(precip, 'precip')
        writeOut(precip, 'tas')
    }
    lapply(models, forModel)
}

mapply(makeDat, names(dirs), dirs, years)