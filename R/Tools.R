"maxNA" <-
function (x) 
{
    return(max(x[!is.na(x)]))
}
"maxNA2" <-
function (x) 
{
    if (is.na(x[1])) 
        return(x[2])
    else {
        if (is.na(x[2])) 
            return(x[1])
        else return(max(x))
    }
}
"meanNA" <-
function (x) 
{
    return(mean(x[!is.na(x)]))
}
"minNA" <-
function (x) 
{
    return(min(x[!is.na(x)]))
}
"minNA2" <-
function (x) 
{
    if (is.na(x[1])) 
        return(x[2])
    else {
        if (is.na(x[2])) 
            return(x[1])
        else return(min(x))
    }
}
#readpartGDAL -------------------------------------------------
#Load geotiff maps partially, only the area of interest 
#delimited by xlim and ylim.
"readpartGDAL" <-
function (x, xlim = NULL, ylim = NULL, ...) 
{
    require(rgdal)
    info <- GDALinfo(x)
    offs <- info[c("ll.x", "ll.y")]
    scl <- info[c("res.x", "res.y")]
    dimn <- info[c("columns", "rows")]
    xs <- seq(offs[1], by = scl[1], length = dimn[1]) + scl[1]/2
    ys <- seq(offs[2], by = scl[2], length = dimn[2]) + scl[2]/2
    xind = 1:length(xs) - 1
    if (!is.null(xlim)) {
        if (!is.numeric(xlim)) 
            stop("xlim must be numeric")
        if (!length(xlim) == 2) 
            stop("xlim must be of length 2")
        if (!diff(xlim) > 0) 
            stop("xlim[1] must be less than xlim[2]")
        xind <- which(xs >= xlim[1] & xs <= xlim[2])
    }
    yind = 1:length(ys)
    if (!is.null(ylim)) {
        if (!is.numeric(ylim)) 
            stop("ylim must be numeric")
        if (!length(ylim) == 2) 
            stop("ylim must be of length 2")
        if (!diff(ylim) > 0) 
            stop("ylim[1] must be less than ylim[2]")
        yind <- which(ys >= ylim[1] & ys <= ylim[2])
    }
    rgdal.offset <- rev(c(min(xind), dimn[2] - max(yind)))
    rgdal.dim <- rev(c(length(xind), length(yind)))
    readGDAL(x, offset = rgdal.offset, region.dim = rgdal.dim, 
        silent = TRUE, ...)
}
"sdNA" <-
function (x) 
{
    return(sd(x[!is.na(x)]))
}

#shapelim -----------------------------------------------------
#Returns the limits of the area of a shape/kml file
"shapelim" <-
function (shapefile, shapedir, shapeext = "shp", around = 0.05) 
{
	while (!tolower(shapeext) %in% c("shp", "kml")) {
        shapeext = readline(cat("Extension is not correct. Please choose between shp and kml. \n"))
        if (shapeext == "") 
            return()
    }
    if (shapeext == "shp") {
        inPoints = readOGR(paste(shapedir, ".", sep = ""), shapefile)
    }
    else {
        inPoints = readOGR(shapedir, shapefile)
    }
    if (dim(coordinates(inPoints))[2] > 2) {
        inPoints = SpatialPointsDataFrame(coords = coordinates(inPoints) [, 1:2], proj4string = CRS(proj4string(inPoints)), data = as.data.frame(inPoints[names(inPoints)]))
    }
    res = c()
    res$xlim = bbox(inPoints)[1, ] + c(-around, around)
    res$ylim = bbox(inPoints)[2, ] + c(-around, around)
    res$proj = CRS(proj4string(inPoints))
    return(res)
}
"timetoMap" <-
function (ndvidirectory, region, year, month, period, type="VITO_CLIP") 
{
    while (!toupper(type) %in% c("GIMMS", "VITO_CLIP", "VITO_VGT")) {
        type = readline(cat("Type is not correct. Please choose between GIMMS, VITO_CLIP and VITO_VGT. \n"))
        if (type == "") 
            stop("Error, type of data NULL.")
    }
    dicperiodVGT = c("01", "11", "21")
    dicmonths = c("01", "02", "03", "04", "05", "06", "07", "08", 
        "09", "10", "11", "12")
    filedate = c("jan15a", "jan15b", "feb15a", "feb15b", "mar15a", 
        "mar15b", "apr15a", "apr15b", "may15a", "may15b", "jun15a", 
        "jun15b", "jul15a", "jul15b", "aug15a", "aug15b", "sep15a", 
        "sep15b", "oct15a", "oct15b", "nov15a", "nov15b", "dec15a", 
        "dec15b")
    codef = c("n07-VIg", "n09-VIg", "n11-VIg", "n14-VIg", "n16-VIg", 
        "n17-VIg")
    if (toupper(type) == "VITO_CLIP") {
        filein = paste(ndvidirectory, region, as.character(year), 
            "M", dicmonths[month], "P", as.character(period - 
                1), ".tif", sep = "")
    }
    if (toupper(type) == "VITO_VGT") {
        filein = paste(ndvidirectory, "NDV_", as.character(year), 
            dicmonths[month], dicperiodVGT[period], "_", region, 
            "_Extract.tif", sep = "")
    }
    if (toupper(type) == "GIMMS") {
        if (period > 2) {
            print("Error gimms data are bimensual.")
        }
        else {
            files = list.files(path = paste(ndvidirectory, ".", 
                sep = ""))
            n = 1
            i = 2 * (month - 1) + period
            filein = paste(ndvidirectory, region, substr(as.character(year), 3, 4), filedate[i], ".", codef[n], ".tif", sep = "")
            while (!(filein %in% paste(ndvidirectory,files,sep=""))) {
                if (n < 7) {
                  n = n + 1
                  filein = paste(ndvidirectory, region, substr(as.character(year), 3, 4), filedate[i], ".", codef[n], ".tif", sep = "")
                }
                else {
                  print(paste("Maps for ", as.character(year), 
                    " - ", filedate[i], " not found.", sep = ""))
                  return()
                }
            }
        }
    }
    if (file.exists(filein)) {
        return(filein)
    }
    else {
        print("Error file doesn't exist !")
        return(filein)
    }
}

#periodtoMap --------------------------------------------------
#Returns the full path of the map of a given period and a given year.
periodtoMap <- 
function (ndvidirectory, region, year, period, type="VITO_CLIP") 
{
    while (!toupper(type) %in% c("GIMMS", "VITO_CLIP", "VITO_VGT")) {
        type = readline(cat("Type is not correct. Please choose between GIMMS, VITO_CLIP and VITO_VGT. \n"))
        if (type == "") 
            stop("Error, type of data NULL.")
    }
    if (toupper(type) == "GIMMS") {
        if (period > 24) {
            stop("Error gimms has only 24 images per year.")
        } else {
        return(timetoMap(ndvidirectory, region, year, ((period-1)%/%2)+1, ((period-1)%%2)+1, type))
        }
    } else {
       if (period > 36) {
            stop("Error vito has only 36 images per year.")
        } else {
        return(timetoMap(ndvidirectory, region, year, ((period-1)%/%3)+1, ((period-1)%%3)+1, type))
        }
    }
}
    
tolist <- 
function(ndvidirectory, region, Ystart, Yend, outfile="list.txt", type = "VITO_CLIP")
{
    period=ifelse(type=="GIMMS",24,36)
    filein=c()
    for (year in Ystart:Yend) {
        for (p in 1:period) {
            filein=c(filein,periodtoMap(ndvidirectory,region,year,p))
        }
    }
    filein=as.data.frame(filein)
    colnames(filein)=length(filein[,1])
    write.table(filein, outfile, quote = FALSE, row.names=FALSE)
}    
