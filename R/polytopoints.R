polytopoints <-
function(shapefile, shapedir, ndvidirectory, region, Ystart, shapeext="shp", type="VITO_CLIP", outname = "buffer", outdir = ".") 
{
    while (!tolower(shapeext) %in% c("shp", "kml")) {
        shapeext = readline(cat("Extension is not correct. Please choose between shp and kml. \n"))
        if (shapeext == "") 
            stop("Error : extension of the shape/kml file is not correct.")
    }
    while (!toupper(type) %in% c("GIMMS", "VITO_CLIP", "VITO_VGT")) {
        type = readline(cat("Type is not correct. Please choose between GIMMS, VITO_CLIP and VITO_VGT. \n"))
        if (type == "") 
            stop("Error : type of data is not correct.")
    }
    if (shapeext == "shp") {
        inPoints = readOGR(paste(shapedir, ".", sep = ""), shapefile)
    }
    else {
        inPoints = readOGR(shapedir, shapefile)
    }
    if (!( length(grep("olygon", class(inPoints))) > 0)) {
        stop("Error : the shape/kml file is not a polygon") 
    }
    if (dim(coordinates(inPoints))[2] > 2) {
        inPoints = SpatialPointsDataFrame(coords = coordinates(inPoints)[, 
            1:2], proj4string = CRS(proj4string(inPoints)), data = as.data.frame(inPoints[names(inPoints)]))
    }
    pro = strsplit(proj4string(inPoints), "[[:punct:]]")[[1]]
    if (!(pro[grep("proj", pro) + 1] == "aea " & pro[grep("ellps", 
        pro) + 1] == "WGS84 ") & type=="GIMMS") {
        inPoints = spTransform(inPoints, CRS("+proj=aea +ellps=WGS84 +lat_1=-19 +lat_2=21 +lat_0=1 +lon_0=20 +x_0=0 +y_0=0"))
    }
    if (!(pro[grep("proj", pro) + 1] == "longlat " & pro[grep("ellps", 
        pro) + 1] == "WGS84 " & pro[grep("datum", pro) + 1] == 
        "WGS84 ") & type!="GIMMS") {
        inPoints = spTransform(inPoints, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
    }
    around=ifelse(type=="GIMMS",16000,0.05)
    filein = timetoMap(ndvidirectory, region, Ystart, 1, 1, type)
    inGrid = readpartGDAL(filein, bbox(inPoints)[1, ] + c(-around, 
        around), bbox(inPoints)[2, ] + c(-around, around))
    a = overlay(inGrid, inPoints)
    xa = coordinates(inGrid)[!is.na(a), 1]
    ya = coordinates(inGrid)[!is.na(a), 2]
    na = c()
    for (i in names(inPoints)) {
        na = cbind(na, as.character(inPoints[[i]][a[!is.na(a)]]))
    }
    colnames(na) = names(inPoints)
    res = SpatialPointsDataFrame(coords = cbind(xa,ya), proj4string = CRS(proj4string(inPoints)), data = as.data.frame(na))
    writeOGR(res, outdir, layer = outname, driver = "ESRI Shapefile")
}
