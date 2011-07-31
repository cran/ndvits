"PeriodAnomaly" <-
function (ndvidirectory, region, Ystart, Yend, period, outname = "anomaly", outext = "show", xlim = NULL, ylim = NULL, type="VITO_CLIP", shapefile = NULL, shapedir = NULL, shapeext = "shp", label = FALSE, pal = "Spectral") 
{
    while (!toupper(type) %in% c("GIMMS", "VITO_CLIP", "VITO_VGT")) {
        type = readline(cat("Type is not correct. Please choose between GIMMS, VITO_CLIP and VITO_VGT. \n"))
        if (type == "") 
            stop("Error : Type is not correct.")
    }
    max=ifelse(type=="GIMMS",10000,255)
    if (!(pal %in% rownames(brewer.pal.info))) {
        stop(paste("Error : pal ", pal , " is not correct name for RColorBrewer palette. \n Please enter 'display.brewer.all()' to see the different choices.",sep="/t"))
    }
    if (!(is.vector(period))) {
        stop("Error : period is not a vector")
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
    while (!tolower(outext) %in% c("jpg", "jpeg", "png", "pdf", "show")) {
        outext = readline(cat("The extension of the output is not correct.\nPlease choose between 'jpg', 'png' or 'pdf' to save the map in the corresponding format or 'show' to see the map directly in R. \n"))
        if (outext == "") 
            return()
    }
    cat("Collecting data from the given period. \n")
    dat=c()
    for (year in seq(Ystart, Yend, 1)) {
        tmp=0
        for (p in period) {
            filein = periodtoMap(ndvidirectory, region, year, p, type) 
            tp = readpartGDAL(filein, xlim, ylim)
            tmp = tmp+tp$band1
        }
        dat = cbind(dat, tmp/(length(period)*max))
        colnames(dat) = c(colnames(dat)[-length(dat[1, ])], as.character(year))
    }
    cat("Computing the mean. \n")
    Mean = apply(dat, 1, meanNA)
    tp = readpartGDAL(timetoMap(ndvidirectory, region, Ystart, 1, 1, type), xlim, ylim)
    tp$band1=Mean*10000
    writeGDAL(tp, paste(outname, "-meanPeriod", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4),".tif", 
            sep = ""), drivername = "GTiff", type = "Int16", 
            mvFlag = -32768)
    Sd = apply(dat, 1, sdNA)
    tp$band1=Sd*10000
    writeGDAL(tp, paste(outname, "-sdPeriod", substr(as.character(Ystart), 3, 4), substr(as.character(Yend), 3, 4),".tif", 
            sep = ""), drivername = "GTiff", type = "Int16", 
            mvFlag = -32768)
    cat("Saving maximum anomaly maps. \n")
    if (outext=="show") {    
        par(ask = TRUE) }
    for (year in seq(Ystart, Yend, 1)) {
        tp$band1 = (dat[, as.character(year)] - Mean) * 10000
        writeGDAL(tp, paste(outname,"P", as.character(year), ".tif", 
            sep = ""), drivername = "GTiff", type = "Int16", 
            mvFlag = -32768)
        savemap(tp, paste(outname, as.character(year), sep = ""), outext, title=paste(outname, " - Period Anomaly ", year, sep=""), shapefile=shapefile, shapedir=shapedir, shapeext=shapeext, label=label, pal=pal)
    }
}
