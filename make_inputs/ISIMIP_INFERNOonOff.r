source("libs/filename.noPath.r")
library(ncdf4)
library(raster)
source("../gitProjectExtras/gitBasedProjects/R/sourceAllLibs.r")
sourceAllLibs("../gitProjectExtras/gitBasedProjects/R/")
sourceAllLibs("../rasterextrafuns/rasterExtras/R/")
source("libs/process_jules_file.r")
source("libs/writeRaster.Standard.r")
options(error=recover)
abcd <- function() source("make_inputs/ISIMIP_INFERNOonOff.r")
countriesMap = raster('../ConFIRE_ISIMIP/data/countries.nc')
ckey = read.csv("../ConFIRE_ISIMIP/data/countries_key.csv")[,2]
genVarID = "genVar-C-cover_clims-sw-precip-ba-tree-fireEmissions2"

fireOffDir = "/hpc/data/d01/hadcam/jules_output/ALL_u-bk886_isimip_0p5deg_origsoil_dailytrif"
fireOnDir = "/hpc/data/d05/cburton/jules_output/u-cf137/"

dirs = list(#historic_on_short  = fireOnDir,
            #historic_off_short = fireOffDir,
            historic_on  = fireOnDir,
            #historic_off = fireOffDir,
            #RCP2.6_on    = fireOnDir,
            RCP6.0_on    = fireOnDir)#,
            #RCP2.6_off   = fireOffDir,
            #RCP6.0_off   = fireOffDir)

years = list(#historic_on_short = 1985:2005,
             #historic_of_short = 1985:2005,
             historic_on  = 1861:2005,
             #historic_off = 1861:2005,
             #RCP2.6_on    = 2006:2099,
             RCP6.0_on    = 2006:2099)#,
             #RCP2.6_off   = 2006:2099,
             #RCP6.0_off   = 2006:2099)

countries = c(#Kenya = 'Kenya', 
              #Indonesia = 'Indonesia',  
              #Malaysia = 'Malaysia', 
              #Brazil = 'Brazil', Paraguay = 'Paraguay', Bolivia = 'Bolivia', 
              #Botswana = 'Botswana', Madagascar = 'Madagascar', 
              #Portugal = 'Portugal', 'NewGuinea' = 'Papua New Guinea', Ghana = 'Ghana', 
              #Russia = 'Russia', Thailand = 'Thailand', IvoryCoast = 'Ivory Coast',
              #Israel = 'Israel', Cambodia = 'Cambodia', Australia = 'Australia', 
              #Canada = 'Canada', USA = 'United States of America', UK = 'United Kingdom',
              Global = NA)

fileIDs = c(temp = "ilamb", sw_down = 'ilamb', npp = "ilamb", smc = "ilamb", 
            cover = "ilamb", cveg = "ilamb", cs_gb = "ilamb", precip = "ilamb", 
            burnt_area = "ilamb", fire_emissions_veg = "ilamb", 
            fire_emissions_dpm = "ilamb", fire_emissions_rpm = "ilamb")

varnames =  c(temp = "t1p5m_gb", sw_down = 'sw_down', npp = "npp_gb", smc = "smc_tot", 
              cover = "frac", cveg = "cv", cs_gb = "cs_gb", precip = "precip", 
              burnt_area = "burnt_area_gb", fire_emissions_veg = "veg_c_fire_emission_gb", 
              fire_emissions_dpm = "burnt_carbon_dpm", fire_emissions_rpm = "burnt_carbon_rpm")

models = c("IPSL-CM5A-LR", "GFDL-ESM2M", "HADGEM2-ES", "MIROC5")[1]

temp_dir = '/data/users/dkelley/ConFIRE_ISIMIP_temp/-makeISIMIPonOffins'
temp_dir_mem = '/data/users/dkelley/ConFIRE_ISIMIP_temp/memSafe/'
out_dir  = '/data/users/dkelley/ConFIRE_ISIMIP/INFERNOonOff/'

coverTypes = list(tallTrees = 1:5, trees = c(1:7), totalVeg = c(1:13), crop = c(10, 12), pas = c(11, 13))
makeDir(out_dir)
try(memSafeFile.remove())
memSafeFile.initialise(temp_dir_mem)
makeDat <- function(id, dir, years_out, out_dir, mask,  extent, country) {
    
    if ( grepl("TS", id) &&  is.null(mask)) return()
    if (!grepl("TS", id) && !is.null(mask)) return()
    print(id)
    print(country)
    years = years_out#c(years_out[1] - 1, years_out)#, tail(years_out, 1) + 1)
    forModel <- function(mod) {
        
        print(id)
        print(mod)
        out_dirM = paste0(out_dir , '/', id)
        makeDir(out_dirM)
        out_dirM = paste(out_dir,  id, mod, '', sep = '/')
        makeDir(out_dirM)
        
        genVarFile = paste0(out_dirM, '/', genVarID, '.Rd')
        if(file.exists(genVarFile)) return()
        tfile0 = paste0(c(temp_dir, country, id, mod, range(years)), collapse = '-')
        dir = paste0(dir, '/', mod, '/')
        files = list.files(dir, full.names = TRUE)
       
        ## select years
        files = files[apply(sapply(years, function(i) grepl(i, files)), 1, any)]
        files = files[substr(files, nchar(files)-2, nchar(files))=='.nc']
       
        openVar <- function(fileID, vname) {
            print(vname)
            tfileC = paste(tfile0 , fileID, vname, '-masked-corrected.Rd', sep = '-')
            print(tfileC)
            
            if (file.exists(tfileC)) {
                load(tfileC)
            } else {
                print(tfileC)
                
                files = files[grepl(fileID, files)]  
                if (substr(id, 1,3) == "RCP") 
                    files = files[grepl(paste0('rcp', substr(id, 4, 6)), files)]
                    #browser()
                
                processSaveFile <- function(file, yr) {
                    dat =  process.jules.file(file, NULL, vname, mask = mask)
                    
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
                
                save(dat, file =  tfileC)
            }
            
            #if (file.exists(tfileC)) {
            #    load(tfileC)
            #} else if(!is.null(mask)) {
            #    countryfy <- function(rs) {                    
            #        applyMask <- function(r) {
            #            r[mask] = NaN
            #            r = raster::crop(r, extent)
            #            r = writeRaster(r, file =memSafeFile())
            #            r
            #       }
            #       if (is.list(rs)) out = lapply(rs, countryfy) else {
            #            fname = filename(rs)
            #            fname = substr(fname, 1, nchar(fname)-3)
            #            fname = paste0(fname, '-', country, '.nc')
            #            print(fname)
            #            if (file.exists(fname)) out = brick(fname)
            #            else {                           
            #                out = layer.apply(rs, applyMask)
            #                writeRaster(out, file = fname, overwrite = TRUE)
            #            }
            #        }
            #        return(out)                  
            #    }
            #    dat = countryfy(dat)
            #    save(dat, file =  tfileC)ø
            #}
            gc()
            return(dat)
        }
        dats = mapply(openVar, fileIDs, varnames)
        
        cover = dats[, 'cover']
        
        #browser()
        #if (mod == "MIROC5") browser()
        makeCover <- function(ty) {
            print(ty)
            strsplitFrac <- function(i) strsplit(filename.noPath(i[[1]], TRUE), 'frac')[[1]][2]
            group <- function(i) {
                ctfile = paste(c(tfile0, 'coverSummed', strsplitFrac(i),
                                 ty, '.nc'), collapse = '-')
                 
                if (file.exists(ctfile)) return(brick(ctfile))
                cv = i[ty]
                
                out = cv[[1]]
                
                for (i in cv[-1]) 
                    out = out + i
                
                                          
                out = writeRaster(out, ctfile, overwrite = TRUE)
                return(out)
            }
            
            ctfile = paste(c(tfile0, 'coverSummed', 
                                 strsplitFrac(cover[[1]]), strsplitFrac(cover[[length(cover)]]),
                                 ty, '.nc'), collapse = '-')
            if (file.exists(ctfile)) return(brick(ctfile))
            
            coverTy = layer.apply(cover, group)
            
            writeRaster(coverTy, ctfile, overwrite = TRUE)
            return(coverTy)
        }
        
        
        covers = lapply(coverTypes, makeCover)
        
        writeOut <- function(dat, name, grab_c = TRUE) {
            file = paste0(out_dirM,  name, '.nc')
            
            if (file.exists(file) && grab_c) return(brick(file))
            print(file)
            dat = dat[[-1]]
            nl = 12*floor(nlayers(dat)/12)
            dat = dat[[1:nl]]
           
            writeRaster.Standard(dat, file)
        }
          
        var2layerWrite <- function(name, file = name) {
            print(name)
            fout = paste0(out_dirM,  file, '.nc')
            
            if (file.exists(fout)) return(brick(fout))
            
            if (length(name) == 1) {
                out = layer.apply(dats[, name], function(i) i)   
            } else {
                forYear <- function(yr) {
                    out = dats[[yr, name[1]]]
                    for (nm in name[-1]) out = out = dats[[yr, nm]]
                    return(out)
                }
                out = layer.apply(1:nrow(dats), forYear)
            }
            
            out = writeOut(out, file)
        }
        
        #fire_emissions_veg = "veg_c_fire_emission_gb",fire_emissions_dpm = "burnt_carbon_dpm", fire_emissions_rpm
        out = list(covers = mapply(writeOut, covers, names(coverTypes), grab_c = FALSE),
                   cveg  = var2layerWrite('cveg'),
                   csoil = var2layerWrite('cs_gb', 'csoil'),
                   temp  = var2layerWrite('temp' ),
                   sw_down  = var2layerWrite('sw_down' ),
                   precip  = var2layerWrite('precip' ),
                   burnt_area  = var2layerWrite('burnt_area' ),
                   fire_emissions = var2layerWrite(c("fire_emissions_veg", 
                                                     "fire_emissions_dpm",
                                                     "fire_emissions_rpm")),
                   npp   = var2layerWrite('npp'  ),
                   smc   = var2layerWrite('smc'  ))
        save(out, file = genVarFile)
        closeAllConnections()
        gc()
    }
    lapply(models, forModel)
    gc()
}

makeCountry <- function(country, cID) {
    if (is.na(cID)) {
        mask = NULL
        extent = NULL
    } else {
        mask = countriesMap!=which(ckey == cID)
        extent = as.vector(apply(xyFromCell(mask, which(mask[] == 0)), 2, range))
        mask = raster::crop(mask,extent)
    }
    out_dirC = paste0(out_dir, '/', country, '/')
    makeDir(out_dirC)
    #extent = as.vector(apply(xyFromCell(mask, mask[] ==1), 2, range))
    mapply(makeDat, names(dirs), dirs, years,
           MoreArgs = list(out_dir = out_dirC, mask = mask, extent = extent,
                           country = country))
}
countriesMap[is.na(countriesMap)] = 0
mapply(makeCountry, names(countries), countries)
memSafeFile.remove()
