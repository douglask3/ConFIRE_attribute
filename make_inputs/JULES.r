
library(ncdf4)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
source("libs/process_jules_file.r")
source("libs/writeRaster.Standard.r")

grab_cache = TRUE

dirs = list(ubl006_WWLLN_FireOff = '/hpc/data/d05/cburton/jules_output/u-bl106_WWLLN_FireOff/')


years = list(historic = 1700:2016)


fileIDs = c(cover = "ilamb", soilM = "Monthly", precip = "ilamb", humid = "ilamb", tas = "ilamb")

varnames =  c(cover = "frac", soilM = "fsmc", precip = "precip", humid = "q1p5m_gb", tas = "t1p5m_gb")

models = c("")

temp_dir = '/data/users/dkelley/ConFIRE_FIREMIP_temp/-makeISIMIPins'
temp_dir_mem = '/data/users/dkelley/ConFIRE_FIREMIP_temp/memSafe/'
out_dir  = '/data/users/dkelley/ConFIRE_FIREMIP/inputs/'
coverTypes = list(trees = c(1:7), totalVeg = c(1:13), crop = c(10, 12), pas = c(11, 13))

memSafeFile.initialise(temp_dir_mem)
makeDat <- function(id, dir, years_out) {
    years = c(years_out[1] - 1, years_out, tail(years_out, 1) + 1)
    forModel <- function(mod) {
        print(id)
        print(mod)
        tfile0 = paste0(c(temp_dir, id, mod, range(years)), collapse = '-')
        dir = paste0(dir, '/', mod, '/')
        
        files = list.files(dir, full.names = TRUE)
       
        ## select years
        files = files[apply(sapply(years, function(i) grepl(i, files)), 1, any)]
        files = files[grepl("365.S3.", files)]
        openVar <- function(fileID, vname) {
            tfile = paste(tfile0 , fileID, vname, '.Rd', sep = '-')
            
            if (file.exists(tfile) && grab_cache) {
                load(tfile)
            } else {
                print(tfile)
                files = files[grepl(fileID, files)]    
                
                processSaveFile <- function(file, yr) {
                    dat =  process.jules.file(file, NULL, vname)
                    if (!is.list(dat)) {
                        dat = writeRaster(dat, paste(tfile0 , fileID, vname, yr,
                                                      '.nc', sep = '-') ,overwrite=TRUE)
                    } else {
                        tfile = paste(tfile0 , fileID, vname, yr, 1:length(dat),
                                      '.nc', sep = '-')
                        dat = mapply(writeRaster, dat, tfile, overwrite = TRUE)
                    }
                    return(dat)
                }
                dat = mapply(processSaveFile, files, years, SIMPLIFY = FALSE)
                
                save(dat, file =  tfile)
            }
            #gc()
            return(dat)
        }
        dats = mapply(openVar, fileIDs, varnames)
        print("yayandwow")

        cover = dats[-1, 'cover']
        #if (mod == "MIROC5") browser()
        makeCover <- function(ty) {
            print(ty)
            group <- function(i) {
                ctfile = paste(c(tfile0, 'coverSummed', ty, '.nc'), collapse = '-')
                if (file.exists(ctfile) && grab_cache) return(brick(ctfile))
                cv = i[ty]
                out = cv[[1]]
                
                for (i in cv[-1]) {
                    print("yay")
                    out = out + i
                }
                               
                out = writeRaster(out, ctfile, overwrite = TRUE)
                return(out)
            }
            
            coverTy = layer.apply(cover, group)
        }
        
        tfile = paste(tfile0, 'Ycovers', '.Rd', sep = '-')
        if (file.exists(tfile) && FALSE) load(tfile)
        else {
            covers = lapply(coverTypes, makeCover)
          
            save(covers, file = tfile)
        }
        #tfile = paste(tfile0, 'soils', '.Rd', sep = '-')
        #if (file.exists(tfile)) load(tfile)
        #else {
            soilM_top    = layer.apply(dats[, 'soilM'], function(i) i[[1]])
            soilM_bottom = layer.apply(dats[-1, 'soilM'], function(i) i[[1]])

            st = 2:(nlayers(soilM_top)-11)
            ed = 13:nlayers(soilM_top)
            soil12 = mapply(function(s, e) sum(soilM_top[[s:e]]), st, ed)
            soil12 = layer.apply(soil12, function(i) i)
            soilM_top =  soilM_top[[-(1:12)]]
            soil12 = soilM_top/soil12
       #     save(soilM_top, soilM_bottom, soil12, file = tfile)
       # }
        precip = layer.apply(dats[-1, 'precip'], function(i) i)
        humid  = layer.apply(dats[-1, 'humid' ], function(i) i)
        tas    = layer.apply(dats[-1, 'tas'   ], function(i) i)

        out_dirM = paste(out_dir,  mod, id, '', sep = '/')
        makeDir(out_dirM)
        writeOut <- function(dat, name) {
            file = paste0(out_dirM,  name, '.nc')
            print(file)
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
        writeOut(humid, 'humid')
        writeOut(tas, 'tas')
        
        gc()
    }
    lapply(models, forModel)
    gc()
}

mapply(makeDat, names(dirs), dirs, years)
memSafeFile.remove()
