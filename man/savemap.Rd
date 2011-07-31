\name{savemap}
\alias{savemap}
\title{
Save a map.
}
\description{
Save a given map in the desire format.
}
\usage{
savemap(map, outname = "map", outext = "show", title = "",
    shapefile = NULL, shapedir = NULL, shapeext = "shp", 
    label = FALSE, pal = "Spectral", ...)
}
\arguments{
  \item{map}{raster to be saved or displayed}
  \item{outname}{name of the file where the map is saved (without extension).}
  \item{outext}{extension of the created map. \cr
If \var{outext} = "show", the map is displayed in the R graphics console.\cr
Choose between "jpg", "png" or "pdf" to save the map in the corresponding format.\cr
By default, \var{outext} = "show".}
  \item{title}{title of the map.\cr
By default, \var{title}="".}
  \item{shapefile}{name of the shape/kml file with objects (points or polygons) to be displayed on the maps.\cr
For a shapefile called \file{name.shp} , \var{shapefile}="name".\cr
For a kml file, write the name of the layer.\cr
If \var{shapefile}=NULL, no object will be added to the map.\cr
By default, \var{shapefile}=NULL.}
  \item{shapedir}{path to the shape/kml file. \cr
For a shapefile, \var{shapedir} is the path to the folder containing the shapefile. 
For example if the directory of the shapefile is \file{C:/Dir/name.shp}, \var{shapefile}="C:/Dir/". \cr
For a kml file, \var{shapedir} is the full directory of the kml file. For example \var{shapefile}="C:/Dir/name.kml". \cr
If \var{shapedir}=NULL, no object will be added to the map. By default, \var{shapedir}=NULL.\cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and "/" in UNIX OS.
}
  \item{shapeext}{extension of the shape/kml file : "shp" for a Esri shapefile or "kml" for a kml file from Google earth .\cr 
By default, \var{shapeext} = "shp".}
  \item{label}{name of the label of the points to be displayed on the map. \cr 
For a shapefile, please indicates the name of the category you want to see.
For a kml file, \var{label}="Name".
If label=FALSE, no label will be plotted. By default, \var{label}=FALSE.}
  \item{pal}{name of the palette of colors used to make the map from the package \pkg{RColorBrewer}.\cr
To see the different palettes available for your map, please type \code{display.brewer.all()}.\cr
By default, \var{pal}="Spectral".}
\item{\dots}{further arguments sent to \command{spplot}. }
}
\author{
Romain Frelat and Bruno Gerard
}
\seealso{
\code{\link{multimap}}
}
\examples{
#local variables
ndvidirectory=paste(system.file("extdata/VITO_Mzimba",
    package="ndvits"), "/", sep="")
region="Mzimba"
#finding the path to the map of the 1st of January 2004
filein=timetoMap(ndvidirectory,region, 2004,1,1)
#loading the raster
map=readpartGDAL(filein)

#executing savemap
savemap(map, outname="Mzimba01012004", title = "NDVI values 1st
    of January 2004 in Mzimba")

#adding points to the map
shape="SLP_Mzimba"
shapedir=paste(system.file("extdata/shape", package="ndvits"),
    "/", sep="")
savemap(map, outname="Mzimba01012004", title = "NDVI values 1st
    of January 2004 in Mzimba", shapefile=shape,
    shapedir=shapedir, label="Village")
}
\keyword{plot}
