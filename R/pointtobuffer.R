"pointtobuffer" <-
function (shapefile, shapedir, ndvidirectory, region, Ystart, shapeext ="shp", outshape = "buffer", outdir = ".", rad = 1, type = "VITO_CLIP") 
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
    coord = c()
    data = c()
    for (i in 1:length(coordinates(inPoints)[, 1])) {
        for (x in -rad:rad) {
            for (y in -rad:rad) {
                coord = rbind(coord, coordinates(inPoints)[i, 
                  ] + c(x * gridparameters(inGrid)[1, 2], y * 
                  gridparameters(inGrid)[2, 2]))
                data = rbind(data, as.data.frame(inPoints[names(inPoints)])[i, 
                  ])
            }
        }
    }
    names(data) = substr(names(data), 1, 10)
    res = SpatialPointsDataFrame(coords = coord, proj4string = CRS(proj4string(inPoints)), 
        data = as.data.frame(data))
    writeOGR(res, outdir, layer = outshape, driver = "ESRI Shapefile")
}
