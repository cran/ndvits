#multimap -----------------------------------------------------
#Display or save multimaps.
"multimap" <-
function (listfiles, names, outname, org, outext = "show", title = "multimap", shapefile = NULL, shapedir = NULL, shapeext = "shp", pal = "Spectral", ...) 
{
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
    if (length(names) != length(listfiles)) {
        cat("Warnings : length of names is different of length of listfiles. \n")
        names = listfiles
    }
    if (org[1]*org[2] < length(listfiles)) {
        cat("Warnings : the layout doesn't allow to show all maps together. \n")
        listfiles=listfiles[1:(org[1]*org[2])]
        names=names[1:(org[1]*org[2])]
    }
    im = readGDAL(listfiles[1])
    dat = im$band1
    for (i in listfiles[-1]) {
        dat = cbind(dat, readGDAL(i)$band1)
    }
    colnames(dat) = paste("X",as.character(1:length(listfiles)),sep="")
    im2 = SpatialGridDataFrame(im, data = as.data.frame(dat))
    sc=mean(apply(bbox(im),1,diff))
    arrow = list("SpatialPolygonsRescale", layout.north.arrow(), offset = c(bbox(im)[1,2]-(sc/5), bbox(im)[2,2]-(sc/5)), scale = sc/7, which = org[1]*org[2], first = FALSE)
    if (typeof(shapefile) == "character" & typeof(shapedir) == "character") {
        if (shapeext == "shp") {
            Bound = readOGR(paste(shapedir, ".", sep = ""), shapefile)
        }
        else {
            Bound = readOGR(shapedir, shapefile)
        }
        if (length(grep("oint", class(Bound))) > 0) {
            lay = list("sp.points", Bound, cex = 0.75, pch = 3, 
                col = "black", first = FALSE)
            spl = list(lay,arrow)
        }
        else {
            if (length(grep("olygon", class(Bound))) > 0) {
                lay = list("sp.polygons", Bound, cex = 0.75, 
                  pch = 3, col = "black", first = FALSE)
                spl = list(lay,arrow)
            }
            else {
                spl = list(arrow)
            }
        }
    }
    else {
        spl = list(arrow)
    }
    if (tolower(outext) == "jpg" || tolower(outext) == "jpeg") {
        jpeg(paste(outname, ".jpg", sep = ""), quality = 90)
        print(spplot(im2, paste("X",as.character(1:length(listfiles)),sep=""), names.attr = names, sp.layout = spl, col.regions = colorRampPalette(brewer.pal(brewer.pal.info[pal,    "maxcolors"], pal)), as.table = TRUE, main = title, layout = org, ...))
        dev.off()
    }
    else {
        if (tolower(outext) == "png") {
            png(paste(outname, ".png", sep = ""))
            print(spplot(im2, paste("X",as.character(1:length(listfiles)),sep=""), names.attr = names, sp.layout = spl, col.regions = colorRampPalette(brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)), as.table = TRUE, main = title, layout = org, ...))
            dev.off()
        }
        else {
            if (tolower(outext) == "pdf") {
                pdf(paste(outname, ".pdf", sep = ""))
                print(spplot(im2, paste("X",as.character(1:length(listfiles)),sep=""), names.attr = names, sp.layout = spl, col.regions = colorRampPalette(brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)), as.table = TRUE, main = title, layout = org, ...))
                dev.off()
            }
            else {
                print(spplot(im2, paste("X",as.character(1:length(listfiles)),sep=""), names.attr = names, sp.layout = spl, col.regions = colorRampPalette(brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)), as.table = TRUE, main = title, layout = org, ...))
            }
        }
    }
}

#savemap ------------------------------------------------------
#Save a given map in the desire format.
"savemap" <-
function (map, outname = "map", outext = "show", title = "", shapefile = NULL, shapedir = NULL, shapeext = "shp", label = FALSE, pal = "Spectral", ...) 
{
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
    sc=mean(apply(bbox(map),1,diff))
    arrow = list("SpatialPolygonsRescale", layout.north.arrow(),offset = c(bbox(map)[1,2]-(sc/10), bbox(map)[2,2]-(sc/10)), scale = sc/15)
    if (typeof(shapefile) == "character" & typeof(shapedir) == "character") {
        if (shapeext == "shp") {
            Bound = readOGR(paste(shapedir, ".", sep = ""), shapefile)
        }
        else {
            Bound = readOGR(shapedir, shapefile)
        }
        if (dim(coordinates(Bound))[2] > 2) {
            Bound = SpatialPointsDataFrame(coords = coordinates(Bound)[, 
                1:2], proj4string = CRS(proj4string(Bound)), 
                data = as.data.frame(Bound[label]))
        }
        if (label != FALSE ) {
        info=names(Bound)
        while (!label %in% info) {
                label = readline(cat(paste("Choose between one of the labels available : \n", list(as.character(info)), "\n : ", sep = "")))
            }
        }
        if (length(grep("oint", class(Bound))) > 0) {
            pt = list("sp.points", Bound, cex = 0.75, pch = 3, 
                col = "black", first = FALSE)
            if (label != FALSE) {
                txt = list("sp.text", coordinates(Bound) - 0.01, 
                  Bound[[label]], first = FALSE)
                lay = list(arrow, pt, txt)
            }
            else {
            lay = list(arrow,pt)
            }
        }
        else {
            if (length(grep("olygon", class(Bound))) > 0) {
                pol = list("sp.polygons", Bound, cex = 0.75, 
                  pch = 3, col = "black", first = FALSE)
                lay = list(arrow,pol)
            }
        }
    }
    else {
        lay = list(arrow)
    }
    if (outext == "jpg" || outext == "jpeg") {
        jpeg(paste(outname, ".jpg", sep = ""), quality = 90)
        print(spplot(map, sp.layout = lay, col.regions = colorRampPalette(brewer.pal(brewer.pal.info[pal, 
            "maxcolors"], pal)), xlab = "Longitude", ylab = "Latitude", 
            add = TRUE, scales = list(draw = TRUE), main = title, 
            font.main = 4))
        dev.off()
    }
    else {
        if (outext == "png") {
            png(paste(outname, ".png", sep = ""))
            print(spplot(map, sp.layout = lay, col.regions = colorRampPalette(brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)), xlab = "Longitude", ylab = "Latitude", add = TRUE, scales = list(draw = TRUE), main = title, font.main = 4))
            dev.off()
        }
        else {
            if (outext == "pdf") {
                pdf(paste(outname, ".pdf", sep = ""))
                print(spplot(map, sp.layout = lay, col.regions = colorRampPalette(brewer.pal(brewer.pal.info[pal, 
                  "maxcolors"], pal)), xlab = "Longitude", 
                  ylab = "Latitude", add = TRUE, scales = list(draw = TRUE), main = title, font.main = 4))
                dev.off()
            }
            else {
                print(spplot(map, sp.layout = lay, col.regions = colorRampPalette(brewer.pal(brewer.pal.info[pal, 
                  "maxcolors"], pal)), xlab = "Longitude", 
                  ylab = "Latitude", add = TRUE, scales = list(draw = TRUE), main = title, font.main = 4))
            }
        }
    }
}
