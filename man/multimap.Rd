\name{multimap}
\alias{multimap}
\title{Create multimap.}
\description{Display or save multimaps.}
\usage{
multimap(listfiles, names, outname, org, outext = "show",
    title = "multimap", shapefile = NULL, shapedir = NULL,
    shapeext = "shp", pal = "Spectral",  ...) 
}
\arguments{
  \item{listfiles}{vector with path to the maps to be plotted together (geotiff files).}
  \item{names}{names of the different maps (display at the top of them).}
  \item{outname}{name of the file where the multi-map is saved (without extension).}
  \item{org}{a two dimensions vector indicating the layout of the maps.\cr
The first number is the number of rows, the second the number of columns.\cr}
  \item{outext}{extension of the created map. \cr
If \var{outext} = "show", the map is displayed in the R graphics console.\cr
Choose between "jpg", "png" or "pdf" to save the map in the corresponding format.\cr
By default, \var{outext} = "show".}
  \item{title}{title of the multi-map.\cr
By default, \var{title}="multimap".}
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
  \item{pal}{name of the palette of colors used to make the map from the package \pkg{RColorBrewer}.\cr
To see the different palettes available for your map, please type \code{display.brewer.all()}.\cr
By default, \var{pal}="Spectral".}
  \item{\dots}{further arguments sent to \command{spplot}. }
}
\author{Romain Frelat and Bruno Gerard.}
\seealso{
\code{\link{savemap}}
}
\examples{
#local variables
ndvidirectory=paste(system.file("extdata/VITO_Mzimba",
    package="ndvits"), "/", sep="")
region="Mzimba"
label=month.abb

#making the list of map to plot together
listname=c()
for (i in 1:12) {
listname=c(listname,timetoMap(ndvidirectory,region, 2004,i,1))
}

#executing multimap function
multimap(listname, label, outname="Mzimba2004", org=c(4,3),
    title = "NDVI values during 2004 in Mzimba")

#adding points to the maps
shape="SLP_Mzimba"
shapedir=paste(system.file("extdata/shape", package="ndvits"),
    "/", sep="")
multimap(listname, label, outname="Mzimba2004", org=c(4,3), 
    title = "NDVI values during 2004 in Mzimba", shapefile=shape,
    shapedir=shapedir)
}
\keyword{plot}
