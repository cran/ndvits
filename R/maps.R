"mapdiff" <-
function (mapin, mapref, mapsd = FALSE, outfile = FALSE) 
{
    if (length(grep("TRUE",bbox(mapin)==bbox(mapref))) != 4) 
        stop("The two maps, mapin and mapref, must have the same size.")
    res = mapref
    if (! is.logical(mapsd)) {
        res$band1 = (mapin$band1 - mapref$band1) / mapsd$band1
    } else {
        res$band1 = mapin$band1 - mapref$band1
    }
    if (outfile != FALSE) {
        writeGDAL(res, outfile, drivername = "GTiff", 
    type = "Int16", mvFlag = -32768)
    }
    return(res)
}
"maplocalstat" <-
function (ndvidirectory, region, Ystart, Yend, xlim = NULL, ylim = NULL, type="VITO_CLIP", outname=FALSE, outext="show", org=c(2,2), shapefile=NULL, shapedir=NULL, shapeext="shp", pal = "Spectral") 
{
    while (!toupper(type) %in% c("GIMMS", "VITO_CLIP", "VITO_VGT")) {
        type = readline(cat("Type is not correct. Please choose between GIMMS, VITO_CLIP and VITO_VGT. \n"))
        if (type == "") 
            stop("Error : Type is not correct.")
    }
    while (!tolower(shapeext) %in% c("shp", "kml")) {
        shapeext = readline(cat("Extension is not correct. Please choose between shp and kml. \n"))
        if (shapeext == "") 
            return()
    }
    if (!(pal %in% rownames(brewer.pal.info))) {
        stop(paste("Error : ", pal , " is not correct name for RColorBrewer palette. \n Please enter 'display.brewer.all()' to see the different choices.", sep=""))
    }
    while (!tolower(outext) %in% c("jpg", "jpeg", "png", "pdf", "show")) {
        outext = readline(cat("The extension of the output is not correct. Please choose between 'jpg', 'png' or 'pdf' to save the map in the corresponding format or 'show' to see the map directly in R. \n"))
        if (outext == "") 
            return()
    }
    if (length(org) != 2) {
        stop("Error : org is not a two-dimensions vector.")
    }
    if (org[1]*org[2] < 4) {
        cat("Error : the layout doesn't display 4 maps together. It is reset to the default value : org=c(2,2).\n")
        org=c(2,2)
    }
    if (!(is.vector(xlim) & length(xlim) == 2 | is.null(xlim))) {
        stop("Error : xlim is not in a convenient format")
    }
    if (!(is.vector(ylim) & length(ylim) == 2 | is.null(ylim))) {
        stop("Error : ylim is not in a convenient format")
    }
    maxNDVI=ifelse(type == "GIMMS",10000,255)
    filein = timetoMap(ndvidirectory, region, Ystart, 1, 1, type)
    inGrid = readpartGDAL(filein, xlim, ylim)
    tp = inGrid
    tp2 = inGrid
    tp3 = inGrid
    tp4 = inGrid
    mat = c()
    for (year in seq(Ystart, Yend, 1)) {
        cat(paste("Processing year ",year,".\n",sep=""))
        for (month in seq(1, 12, 1)) {
            for (period in seq(1, 3, 1)) {
                if (!(type=="GIMMS" & period==3)) {
                    filein = timetoMap(ndvidirectory, region,  year, month, period, type)
                    inGrid = readpartGDAL(filein, xlim, ylim)
                    mat = cbind(mat, inGrid$band1/maxNDVI)
                }
            }
        }
    }
    tp$band1 = apply(mat, 1, meanNA)*10000
    tp2$band1 = apply(mat, 1, maxNA)*10000
    tp3$band1 = apply(mat, 1, minNA)*10000
    tp4$band1 = apply(mat, 1, sdNA)*10000
    res = c()
    res$mean = tp
    res$max = tp2
    res$min = tp3
    res$sd = tp4
    if (outname != FALSE) {
        writeGDAL(tp, paste(outname, "mean", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4), ".tif", sep = ""), drivername = "GTiff", type = "Int16", mvFlag = -32768)
        writeGDAL(tp2, paste(outname, "max", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4), ".tif", sep = ""), drivername = "GTiff", type = "Int16", mvFlag = -32768)
        writeGDAL(tp3, paste(outname, "min", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4), ".tif", sep = ""), drivername = "GTiff", type = "Int16", mvFlag = -32768)
        writeGDAL(tp4, paste(outname, "sd", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4), ".tif", sep = ""), drivername = "GTiff", type = "Int16", mvFlag = -32768)
        #creating the multimap
        listfile=paste(outname, c("max","min","mean","sd"), substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4), ".tif", sep = "")
        multimap(listfile, c("max","min","mean","sd"), outname=paste(outname, "Stat", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4), sep = ""), org=org, outext=outext, title=paste(outname, " NDVI local statistics ",as.character(Ystart), " - ", as.character(Yend), sep=""), shapefile=shapefile, shapedir=shapedir, shapeext=shapeext, pal=pal)
    }
    return(res)
}

"mapmaxyear" <-
function (ndvidirectory, region, year, xlim = NULL, ylim = NULL, 
    outfile = FALSE, type="VITO_CLIP") 
{
    while (!toupper(type) %in% c("GIMMS", "VITO_CLIP", "VITO_VGT")) {
        type = readline(cat("Type is not correct. Please choose between GIMMS, VITO_CLIP and VITO_VGT. \n"))
        if (type == "") 
            stop("Error : Type is not correct.")
    }
    if (!(is.vector(xlim) & length(xlim) == 2 | is.null(xlim))) {
        stop("Error : xlim is not in a convenient format")
    }
    if (!(is.vector(ylim) & length(ylim) == 2 | is.null(ylim))) {
        stop("Error : ylim is not in a convenient format")
    }
    maxNDVI=ifelse(type == "GIMMS",10000,255)
    filein = timetoMap(ndvidirectory, region, year, 1, 1, type)
    inGrid = readpartGDAL(filein, xlim, ylim)
    tp = inGrid
    for (month in seq(1, 12, 1)) {
        for (period in seq(1, 3, 1)) {
            if (!(type=="GIMMS" & period==3)) {
                filein = timetoMap(ndvidirectory, region,  year, month, period, type)
                inGrid = readpartGDAL(filein, xlim, ylim)
                tp$band1 = apply(cbind(tp$band1, inGrid$band1), 1, 
                maxNA2)
            }
        }
    }
    tp$band1=(tp$band1/maxNDVI)*10000
    if (outfile != FALSE) {
        writeGDAL(tp, outfile, drivername = "GTiff", type = "Int16", 
        mvFlag = -32768)
    }
    return(tp)
}
