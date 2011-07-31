\name{shapelim}
\alias{shapelim}
\title{
limits of a shape/kml file
}
\description{
Returns the limits of the area of a shape/kml file}
\usage{
shapelim(shapefile, shapedir, shapeext = "shp", around = 0.05)
}
\arguments{
  \item{shapefile}{name of the shape/kml file. \cr
For a shapefile called \file{name.shp} , \var{shapefile}="name".\cr
For a kml file, write the name of the layer.}
  \item{shapedir}{path to the shape/kml file. \cr
For a shapefile, \var{shapedir} is the path to the folder containing the shapefile. 
For example if the directory of the shapefile is \file{C:/Dir/name.shp}, \var{shapefile}="C:/Dir/". \cr
For a kml file, \var{shapedir} is the full directory of the kml file. For example \var{shapefile}="C:/Dir/name.kml". \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and "/" in UNIX OS.}
  \item{shapeext}{
extension of the shape/kml file : "shp" for a Esri shapefile or "kml" for a kml file from Google earth .\cr 
By default, \var{shapeext} = "shp".
}
  \item{around}{
extra space added around the area of the shapefile. By default, \var{around}=0.05.\cr
Please note that \var{around} have to be in the same units and projection as the ndvi maps.
}
}
\value{
Returns a list with :
\item{xlim}{a two dimensions vector with the longitude borders of the shapefile}
\item{ylim}{a two dimensions vector with the latitude borders of the shapefile}
\item{proj}{string containing the projection parameters of the shapefile}
}
\author{
Romain Frelat and Bruno Gerard
}
\examples{
#with a shape file :
shape="SLP_Mzimba"
shapedir=paste(system.file("extdata/shape", package="ndvits"), "/", sep="")
lim1=shapelim(shape,shapedir)
print(lim1)

#with a kml file
shape="SLP-SAs_CropLand"
shapedir=system.file("extdata/shape/SLP-SAs.kml", package="ndvits")
lim2=shapelim(shape,shapedir,shapeext="kml")
print(lim2)
}
\keyword{spatial}
