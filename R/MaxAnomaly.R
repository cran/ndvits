#MaxAnomaly ---------------------------------------------------
#
"MaxAnomaly" <-
function (ndvidirectory, region, Ystart, Yend, outname = "anomaly", outext = "show", xlim = NULL, ylim = NULL, type="VITO_CLIP", shapefile = NULL, shapedir = NULL, shapeext = "shp", label = FALSE, pal = "Spectral") 
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
    while (!tolower(shapeext) %in% c("shp", "kml")) {
        shapeext = readline(cat("Extension is not correct. Please choose between shp and kml. \n"))
        if (shapeext == "") 
            return()
    }
    if (!(pal %in% rownames(brewer.pal.info))) {
        stop(paste("Error : ", pal , " is not correct name for RColorBrewer palette. \n Please enter 'display.brewer.all()' to see the different choices.", sep=""))
    }
    while (!tolower(outext) %in% c("jpg", "jpeg", "png", "pdf", "show")) {
        outext = readline(cat("The extension of the output is not correct.\nPlease choose between 'jpg', 'png' or 'pdf' to save the map in the corresponding format or 'show' to see the map directly in R. \n"))
        if (outext == "") 
            return()
    }
    dat = c()
    tab = c()
    cat("Computing the maximum of each year. \n")
    for (year in seq(Ystart, Yend, 1)) {
        cat(paste("Processing year ",year,".\n",sep=""))
        dat = cbind(dat, mapmaxyear(ndvidirectory, region, 
            year, xlim, ylim, type=type)$band1)
        colnames(dat) = c(colnames(dat)[-length(dat[1, ])], as.character(year))
    }
    cat("Computing the mean. \n")
    Mean = apply(dat, 1, meanNA)
    tp = readpartGDAL(timetoMap(ndvidirectory, region, Ystart, 1, 1, type), xlim, ylim)
    tp$band1=Mean
    writeGDAL(tp, paste(outname, "-meanMax", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4),".tif", 
            sep = ""), drivername = "GTiff", type = "Int16", 
            mvFlag = -32768)
    Sd = apply(dat, 1, sdNA)
    tp$band1=Sd
    writeGDAL(tp, paste(outname, "-sdMax", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4),".tif", 
            sep = ""), drivername = "GTiff", type = "Int16", 
            mvFlag = -32768)
    cat("Saving maximum anomaly maps. \n")
    if (outext=="show") {    
        par(ask = TRUE) }
    for (year in seq(Ystart, Yend, 1)) {
        tp$band1 = (dat[, as.character(year)] - Mean)
        writeGDAL(tp, paste(outname,"M", as.character(year), ".tif", 
            sep = ""), drivername = "GTiff", type = "Int16", 
            mvFlag = -32768)
        savemap(tp, paste(outname, as.character(year), sep = ""), outext, title=paste(outname, " - Maximum Anomaly ", year, sep=""), shapefile=shapefile, shapedir=shapedir, shapeext=shapeext, label=label, pal=pal)
    }
}
