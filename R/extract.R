"ExtractFile" <-
function (shapefile, shapedir, listfile, outfile, period, shapeext = "shp") 
{
    list=read.table(listfile, header = FALSE, sep="\t")
    nts=as.numeric(as.character(list[1,1]))
    if (is.na(nts) | length(list[,1])<nts+1) {
        stop("The file is not in the good format. The first line should contain the number of images to be processed.")
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
    pro = strsplit(proj4string(inPoints), "[[:punct:]]")[[1]]
    inGrid=readGDAL(list[2,1])
    proI = strsplit(proj4string(inGrid), "[[:punct:]]")[[1]]
    if (!(pro[grep("proj", pro) + 1] == proI[grep("proj", proI) + 1] & pro[grep("ellps", pro) + 1] == proI[grep("ellps", proI) + 1])) {
        inPoints = spTransform(inPoints, CRS(proj4string(inGrid)))
    }
    ndvits=c()
    for (filein in list[2:(nts+1),1]) {
        inGrid = readpartGDAL(filein, bbox(inPoints)[1,] + c(-gridparameters(inGrid)[1, 2], gridparameters(inGrid)[1, 2]), bbox(inPoints)[2, ] + c(-gridparameters(inGrid)[2, 2], gridparameters(inGrid)[2, 2]))
        if (length(grep("oint", class(inPoints))) > 0) {
            ndvits = cbind(ndvits, overlay(inGrid, inPoints)$band1)
        } else {
            if (length(grep("olygon", class(inPoints))) > 0) {
                ndvits = cbind(ndvits, inGrid@data[!is.na(overlay(inGrid, inPoints)), ])
            } else {
                print("The type of ths shp/kml file is unknown.")
            }   
        }
    }
    if (shapeext == "shp") {
        if (length(grep("oint", class(inPoints))) > 0) {
            all = data.frame(inPoints[names(inPoints)], ndvits)
        }
        if (length(grep("olygon", class(inPoints))) > 0) {
            a = overlay(inGrid, inPoints)
            xa = coordinates(inGrid)[!is.na(a), 1]
            ya = coordinates(inGrid)[!is.na(a), 2]
            na = c()
            for (i in names(inPoints)) {
                na = cbind(na, as.character(inPoints[[i]][a[!is.na(a)]]))
            }
            colnames(na) = names(inPoints)
            all = data.frame(cbind(xa, ya, na, ndvits))
        }
    }
    else {
        if (length(grep("oint", class(inPoints))) > 0) {
            all = data.frame(cbind(as.character(inPoints$Name), 
                ndvits))
        }
        if (length(grep("olygon", class(inPoints))) > 0) {
            a = overlay(inGrid, inPoints)
            xa = coordinates(inGrid)[!is.na(a), 1]
            ya = coordinates(inGrid)[!is.na(a), 2]
            name = as.character(inPoints$Name[a[!is.na(a)]])
            all = data.frame(cbind(name, ndvits))
        }
    }
    nyear=nts%/%period
    if (nts%%period !=0) {
        cat("Warning message: \nThe time series doesn't stop with at the end of a year, the last year is not complete.\n") 
    }
    write(paste(as.character(nyear), as.character(period), as.character(nts),sep=" "), outfile, sep = "\t")
    write.table(all, outfile, append = TRUE, quote = FALSE, row.names = TRUE, sep = "\t")
    return(all)
}

"ExtractGIMMS" <-
function (shapefile, shapedir, ndvidirectory, region, outfile, 
    Ystart, Yend, shapeext = "shp") 
{
    while (!toupper(region) %in% c("AF", "AZ", "EA", "NA", "SA", 
        "")) {
        region = readline(cat("Region is not correct. Please choose between AF, AZ, EA, NA and SA. \n"))
    }
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
        inPoints = SpatialPointsDataFrame(coords = coordinates(inPoints)[, 
            1:2], proj4string = CRS(proj4string(inPoints)), data = as.data.frame(inPoints[names(inPoints)]))
    }
    pro = strsplit(proj4string(inPoints), "[[:punct:]]")[[1]]
    if (!(pro[grep("proj", pro) + 1] == "aea " & pro[grep("ellps", 
        pro) + 1] == "WGS84 ")) {
        inPoints = spTransform(inPoints, CRS("+proj=aea +ellps=WGS84 +lat_1=-19 +lat_2=21 +lat_0=1 +lon_0=20 +x_0=0 +y_0=0"))
    }
    filedate = c("jan15a", "jan15b", "feb15a", "feb15b", "mar15a", 
        "mar15b", "apr15a", "apr15b", "may15a", "may15b", "jun15a", 
        "jun15b", "jul15a", "jul15b", "aug15a", "aug15b", "sep15a", 
        "sep15b", "oct15a", "oct15b", "nov15a", "nov15b", "dec15a", 
        "dec15b")
    codef = c("n07-VIg", "n09-VIg", "n11-VIg", "n14-VIg", "n16-VIg", 
        "n17-VIg")
    ndvits = c()
    files = list.files(path = paste(ndvidirectory, ".", sep = ""))
    for (year in seq(Ystart, Yend, 1)) {
        cat(paste("Processing year", as.character(year), "\n"))
        for (i in seq(1, 24, 1)) {
            n = 1
            filein = paste(region, substr(as.character(year), 
                3, 4), filedate[i], ".", codef[n], ".tif", sep = "")
            while (!(filein %in% files)) {
                if (n < 7) {
                  n = n + 1
                  filein = paste(region, substr(as.character(year), 
                    3, 4), filedate[i], ".", codef[n], ".tif", 
                    sep = "")
                }
                else {
                  print(paste("Maps for ", as.character(year), 
                    " - ", filedate[i], " not found.", sep = ""))
                  return()
                }
            }
            inGrid = readpartGDAL(paste(ndvidirectory, filein, 
                sep = ""), bbox(inPoints)[1, ] + c(-25000, 25000), 
                bbox(inPoints)[2, ] + c(-25000, 25000))
            if (length(grep("oint", class(inPoints))) > 0) {
                ndvits = cbind(ndvits, overlay(inGrid, inPoints)$band1)
                colnames(ndvits) = c(colnames(ndvits)[-length(ndvits[1, 
                  ])], paste(as.character(year), "_", as.character(i), 
                  sep = ""))
            }
            else {
                if (length(grep("olygon", class(inPoints))) > 
                  0) {
                  ndvits = cbind(ndvits, inGrid@data[!is.na(overlay(inGrid, 
                    inPoints)), ])
                  colnames(ndvits) = c(colnames(ndvits)[-length(ndvits[1, 
                    ])], paste(as.character(year), "_", as.character(i), 
                    sep = ""))
                }
                else {
                  print(paste("The type of the ", shapeext, " file is unknown.", 
                    sep = ""))
                }
            }
        }
    }
    if (shapeext == "shp") {
        if (length(grep("oint", class(inPoints))) > 0) {
            all = data.frame(inPoints[names(inPoints)], ndvits)
        }
        if (length(grep("olygon", class(inPoints))) > 0) {
            a = overlay(inGrid, inPoints)
            xa = coordinates(inGrid)[!is.na(a), 1]
            ya = coordinates(inGrid)[!is.na(a), 2]
            na = c()
            for (i in names(inPoints)) {
                na = cbind(na, as.character(inPoints[[i]][a[!is.na(a)]]))
            }
            colnames(na) = names(inPoints)
            all = data.frame(cbind(xa, ya, na, ndvits))
        }
    }
    else {
        if (length(grep("oint", class(inPoints))) > 0) {
            all = data.frame(cbind(as.character(inPoints$Name), 
                ndvits))
        }
        if (length(grep("olygon", class(inPoints))) > 0) {
            a = overlay(inGrid, inPoints)
            xa = coordinates(inGrid)[!is.na(a), 1]
            ya = coordinates(inGrid)[!is.na(a), 2]
            name = as.character(inPoints$Name[a[!is.na(a)]])
            all = data.frame(cbind(name, ndvits))
        }
    }
    write.table(all, outfile, quote = FALSE, row.names = TRUE, sep = "\t")
    return(all)
}
"ExtractVGT" <-
function (shapefile, shapedir, ndvidirectory, region, outfile = "TS.txt", 
    Ystart, Yend, shapeext = "shp") 
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
        inPoints = SpatialPointsDataFrame(coords = coordinates(inPoints)[, 
            1:2], proj4string = CRS(proj4string(inPoints)), data = as.data.frame(inPoints[names(inPoints)]))
    }
    pro = strsplit(proj4string(inPoints), "[[:punct:]]")[[1]]
    if (!(pro[grep("proj", pro) + 1] == "longlat " & pro[grep("ellps", 
        pro) + 1] == "WGS84 " & pro[grep("datum", pro) + 1] == 
        "WGS84 ")) {
        inPoints = spTransform(inPoints, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
    }
    ndvits = c()
    dicperiod = c("01", "11", "21")
    dicmonths = c("01", "02", "03", "04", "05", "06", "07", "08", 
        "09", "10", "11", "12")
    for (year in seq(Ystart, Yend, 1)) {
        cat(paste("Processing year", as.character(year), "\n"))
        for (month in seq(1, 12, 1)) {
            for (period in seq(1, 3, 1)) {
                filein = paste(ndvidirectory, "NDV_", as.character(year), 
                  dicmonths[month], dicperiod[period], "_", region, 
                  "_Extract.tif", sep = "")
                inGrid = readpartGDAL(filein, bbox(inPoints)[1, 
                  ] + c(-0.05, 0.05), bbox(inPoints)[2, ] + c(-0.05, 
                  0.05))
                if (length(grep("oint", class(inPoints))) > 0) {
                  ndvits = cbind(ndvits, overlay(inGrid, inPoints)$band1)
                  colnames(ndvits) = c(colnames(ndvits)[-length(ndvits[1, 
                    ])], paste(as.character(year), dicmonths[month], 
                    as.character(period), sep = ""))
                }
                else {
                  if (length(grep("olygon", class(inPoints))) > 
                    0) {
                    ndvits = cbind(ndvits, inGrid@data[!is.na(overlay(inGrid, 
                      inPoints)), ])
                    colnames(ndvits) = c(colnames(ndvits)[-length(ndvits[1, 
                      ])], paste(as.character(year), dicmonths[month], 
                      as.character(period), sep = ""))
                  }
                  else {
                    print("The type of ths shp/kml file is unknown.")
                  }
                }
            }
        }
    }
    if (shapeext == "shp") {
        if (length(grep("oint", class(inPoints))) > 0) {
            all = data.frame(inPoints[names(inPoints)], ndvits)
        }
        if (length(grep("olygon", class(inPoints))) > 0) {
            a = overlay(inGrid, inPoints)
            xa = coordinates(inGrid)[!is.na(a), 1]
            ya = coordinates(inGrid)[!is.na(a), 2]
            na = c()
            for (i in names(inPoints)) {
                na = cbind(na, as.character(inPoints[[i]][a[!is.na(a)]]))
            }
            colnames(na) = names(inPoints)
            all = data.frame(cbind(xa, ya, na, ndvits))
        }
    }
    else {
        if (length(grep("oint", class(inPoints))) > 0) {
            all = data.frame(cbind(as.character(inPoints$Name), 
                ndvits))
        }
        if (length(grep("olygon", class(inPoints))) > 0) {
            a = overlay(inGrid, inPoints)
            xa = coordinates(inGrid)[!is.na(a), 1]
            ya = coordinates(inGrid)[!is.na(a), 2]
            name = as.character(inPoints$Name[a[!is.na(a)]])
            all = data.frame(cbind(name, ndvits))
        }
    }
    write.table(all, outfile, quote = FALSE, row.names = TRUE, sep = "\t")
    return(all)
}
"ExtractVito" <-
function (shapefile, shapedir, ndvidirectory, region, outfile = "TS.txt", 
    Ystart, Yend, shapeext = "shp") 
{
    while (!tolower(shapeext) %in% c("shp", "kml")) {
        shapeext = readline(cat("Extension is not correct. Please choose between shp and kml. \n"))
        if (shapeext == "") 
            return()
    }
    if (tolower(shapeext) == "shp") {
        inPoints = readOGR(paste(shapedir, ".", sep = ""), shapefile)
    }
    else {
        inPoints = readOGR(shapedir, shapefile)
        name = inPoints$Name
    }
    if (dim(coordinates(inPoints))[2] > 2) {
        inPoints = SpatialPointsDataFrame(coords = coordinates(inPoints)[, 
            1:2], proj4string = CRS(proj4string(inPoints)), data = as.data.frame(inPoints[names(inPoints)]))
    }
    pro = strsplit(proj4string(inPoints), "[[:punct:]]")[[1]]
    if (!(pro[grep("proj", pro) + 1] == "longlat " & pro[grep("ellps", 
        pro) + 1] == "WGS84 " & pro[grep("datum", pro) + 1] == 
        "WGS84 ")) {
        inPoints = spTransform(inPoints, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84"))
    }
    ndvits = c()
    for (year in seq(Ystart, Yend, 1)) {
        cat(paste("Processing year", as.character(year), "\n"))
        for (month in seq(1, 12, 1)) {
            if (month < 10) {
                monthstring = paste("0", as.character(month), 
                  sep = "")
            }
            else {
                monthstring = as.character(month)
            }
            for (period in seq(0, 2, 1)) {
                filein = paste(ndvidirectory, region, as.character(year), 
                  "M", monthstring, "P", as.character(period), 
                  ".tif", sep = "")
                inGrid = readpartGDAL(filein, bbox(inPoints)[1, 
                  ] + c(-0.05, 0.05), bbox(inPoints)[2, ] + c(-0.05, 
                  0.05))
                if (length(grep("oint", class(inPoints))) > 0) {
                  ndvits = cbind(ndvits, overlay(inGrid, inPoints)$band1)
                  colnames(ndvits) = c(colnames(ndvits)[-length(ndvits[1, 
                    ])], paste(as.character(year), monthstring, 
                    as.character(period), sep = ""))
                }
                else {
                  if (length(grep("olygon", class(inPoints))) > 
                    0) {
                    ndvits = cbind(ndvits, inGrid@data[!is.na(overlay(inGrid, 
                      inPoints)), ])
                    colnames(ndvits) = c(colnames(ndvits)[-length(ndvits[1, 
                      ])], paste(as.character(year), monthstring, 
                      as.character(period), sep = ""))
                  }
                  else {
                    print("The type of ths shp/kml file is unknown.")
                  }
                }
            }
        }
    }
    if (shapeext == "shp") {
        if (length(grep("oint", class(inPoints))) > 0) {
            all = data.frame(inPoints[names(inPoints)], ndvits)
        }
        if (length(grep("olygon", class(inPoints))) > 0) {
            a = overlay(inGrid, inPoints)
            xa = coordinates(inGrid)[!is.na(a), 1]
            ya = coordinates(inGrid)[!is.na(a), 2]
            na = c()
            for (i in names(inPoints)) {
                na = cbind(na, as.character(inPoints[[i]][a[!is.na(a)]]))
            }
            colnames(na) = names(inPoints)
            all = data.frame(cbind(xa, ya, na, ndvits))
        }
    }
    else {
        if (length(grep("oint", class(inPoints))) > 0) {
            all = data.frame(cbind(as.character(inPoints$Name), 
                ndvits))
        }
        if (length(grep("olygon", class(inPoints))) > 0) {
            a = overlay(inGrid, inPoints)
            xa = coordinates(inGrid)[!is.na(a), 1]
            ya = coordinates(inGrid)[!is.na(a), 2]
            name = as.character(inPoints$Name[a[!is.na(a)]])
            all = data.frame(cbind(name, ndvits))
        }
    }
    write.table(all, outfile, quote = FALSE, row.names = TRUE, sep = "\t")
    return(all)
}
"normNDVI" <-
function (TS, maxNDVI) 
{
    ndvi = TS[grep("X", names(TS))]
    for (i in 1:length(ndvi)) {
        ndvi[, i] = as.numeric(as.character(ndvi[, i]))/maxNDVI
    }
    return(ndvi)
}
"STLperArea" <-
function (TS, area, outfile = FALSE, Ystart, period = 36, fct = "mean", 
    SGfilter = TRUE, nSG = "5,5", DSG = 0) 
{
    if (!require(RTisean)) {
        return()
    }
    while (!tolower(fct) %in% c("mean", "max", "min")) {
        fct = readline(cat("Function is not correct. Please choose between mean, max and min."))
        if (fct == "") 
            return()
    }
    finalndvi = c()
    area = as.factor(area)
    if (outfile != FALSE) {
        pdf(outfile)
    }
    else {
        par(ask = TRUE)
    }
    for (i in unique(area)) {
        par(oma = c(0, 0, 3, 0))
        ndviLoc = TS[area == i, ]
        ndviM = apply(ndviLoc, 2, fct)
        plot(ts(ndviM, start = Ystart, freq = period), ylim = c(0, 
            1), type = "l", xlab = "time", ylab = "NDVI", main = 
            paste("NDVI Time Serie -", fct, "over", length(ndviLoc[, 1]), "points"))
        mtext(side = 3, text = i, outer = TRUE, cex = 1.6)
        if (SGfilter) {
            filtndvi = sav_gol(ndviM, n = nSG, D = DSG)
            ndviMSG = ts(filtndvi[, 1], start = Ystart, freq = period)
            lines(ndviMSG, col = "red")
            if (min(ndviMSG) > 0.2 | max(ndviMSG) > 0.85) {
                legend("bottomright", legend = c(paste(fct, "NDVI", 
                  sep = " "), "SG filter"), col = c("black", 
                  "red"), lwd = 2)
            }
            else {
                legend("topright", legend = c(paste(fct, "NDVI", 
                  sep = " "), "SG filter"), col = c("black", 
                  "red"), lwd = 2)
            }
            filt = "SG Filtered"
        }
        else {
            filt = ""
            ndviMSG = ts(ndviM, start = Ystart, freq = period)
        }
        title = c(i, paste(fct, filt, "NDVI time series"))
        plot(stl(ndviMSG, s.window = "periodic"), main = title, robust=TRUE, na.action=na.approx)
        finalndvi = rbind(finalndvi, ndviMSG)
        rownames(finalndvi) = c(rownames(finalndvi)[-length(finalndvi[, 
            1])], i)
    }
    if (outfile != FALSE) {
        dev.off()
    }
    colnames(finalndvi) = as.character(round(time(ndviMSG), 2))
    return(finalndvi)
}

"TimeSeriesAnalysis" <-
function (shapefile, shapedir, ndvidirectory, region, Ystart, 
    Yend, outfile = "TS.txt", outfile2 = "TS.pdf", outfile3 = FALSE, 
    shapeext = "shp", fct = "mean", SGfilter = TRUE, nSG = "5,5", 
    DSG = 0, title = "NDVI time series", type = "VITO_CLIP", nb = NULL) 
{
    while (!tolower(shapeext) %in% c("shp", "kml")) {
        shapeext = readline(cat("Extension is not correct. Please choose between shp and kml. \n"))
        if (shapeext == "") 
            return()
    }
    if (tolower(shapeext) == "shp") {
        inPoints = readOGR(paste(shapedir, ".", sep = ""), shapefile)
        info = names(inPoints)
        if (length(info) == 1) {
            fac = as.factor(inPoints[[info]])
        }
        else {
            f = ""
            while (!f %in% info) {
                f = readline(cat(paste("choose between one of the grouping factor available : \n", 
                  list(as.character(info)), "\n : ", sep = "")))
            }
            fac = as.factor(inPoints[[f]])
        }
    }
    else {
        inPoints = readOGR(shapedir, shapefile)
        fac = as.factor(inPoints$Name)
        f="name"
    }
    while (!toupper(type) %in% c("GIMMS", "VITO_CLIP", "VITO_VGT", 
        "TEXT", "FILES")) {
        type = readline(cat("Type is not correct. Please choose between GIMMS, VITO_CLIP, VITO_VGT, TEXT and FILES. \n"))
        if (type == "") 
            return()
    }
    if (toupper(type) == "VITO_CLIP") {
        TS = ExtractVito(shapefile, shapedir, ndvidirectory, 
            region, outfile, Ystart, Yend, shapeext)
        period = 36
        max = 255
    }
    if (toupper(type) == "VITO_VGT") {
        TS = ExtractVGT(shapefile, shapedir, ndvidirectory, region, 
            outfile, Ystart, Yend, shapeext)
        period = 36
        max = 255
    }
    if (toupper(type) == "GIMMS") {
        TS = ExtractGIMMS(shapefile, shapedir, ndvidirectory, 
            region, outfile, Ystart, Yend, shapeext)
        period = 24
        max = 10000
    }
    if (toupper(type) == "TEXT") {
        period = as.numeric(readline(cat("how many observation per year ?\n")))
        max = as.numeric(readline(cat("value maximum ?\n")))
        TS = read.table(ndvidirectory, header = TRUE, sep = "\t", 
            row.names = 1)
    }
    if (toupper(type) == "FILES") {
        period = as.numeric(readline(cat("how many observation per year ?\n")))
        max = as.numeric(readline(cat("value maximum ?\n")))
        TS = ExtractFile(shapefile, shapedir, ndvidirectory, outfile, period, shapeext)
    }
    if (length(grep("olygon", class(inPoints))) > 0) {
        while (!f %in% names(TS)) {
                f = readline(cat(paste("error choose between : \n", 
                  list(as.character(head(names(TS)))), "\n : ", sep = "")))
        }
        fac = as.factor(TS[,f])
    }
    ndvi = normNDVI(TS, max)
    TS2 = STLperArea(ndvi, fac, outfile2, Ystart, 
        period, fct, SGfilter, nSG, DSG)
    PlotTS(TS2, outfile3, Ystart, period, title, nb)
    return(TS2)
}
