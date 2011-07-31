\name{ExtractFile}
\alias{ExtractFile}
\title{
Extract time series from a list of images
}
\description{
Extract time series from a file containing a list of images
}
\usage{
ExtractFile(shapefile, shapedir, listfile, outfile, period, 
    shapeext = "shp")
}
\arguments{
  \item{shapefile}{name of the shape/kml file with multi-points or polygons from where time series are extracted. \cr
For a shapefile called \file{name.shp} , \var{shapefile}="name".\cr
For a kml file, write the name of the layer.}
  \item{shapedir}{path to the shape/kml file. \cr
For a shapefile, \var{shapedir} is the path to the folder containing the shapefile. 
For example if the directory of the shapefile is \file{C:/Dir/name.shp}, \var{shapefile}="C:/Dir/". \cr
For a kml file, \var{shapedir} is the full directory of the kml file. For example \var{shapefile}="C:/Dir/name.kml". \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and "/" in UNIX OS.}
  \item{listfile}{file containing a list of images}
  \item{outfile}{name of the file where the time series are saved (extension ".txt").\cr By default, \var{outfile} = " .".}
  \item{period}{number of observation per year}
  \item{shapeext}{extension of the shape/kml file : "shp" for a Esri shapefile or "kml" for a kml file from Google earth .\cr 
By default, \var{shapeext} = "shp".}
}
\details{
\var{listfile} should respect the following format :\cr
\verb{N}\cr
\verb{path\image1.tif}\cr
\verb{path\image2.tif}\cr
\verb{...}\cr
\verb{path\imageN.tif}\cr
with \verb{N} the number of image files. \cr

The list of images should be in chronological order, and include the complete path to the file, unless the files are located in the working directory from which R is run. Only the N first images in the file will be processed.
}
\value{Return the time series of the points/polygons in the shapefile over images in the shapefile.}
\author{Romain Frelat and Bruno Gerard}
\section{Warning }{It is the responsability of the user to have complete and ordered dataset : list of images have to be in chronological order and without missing image.}
\seealso{\code{\link{TimeSeriesAnalysis}},\code{\link{tolist}}}
\examples{
#Be carefull : it is only an example !
#It is not the good way to extract time series
#when you have the maps in an accepted nomenclature.
#You should always prefer ExtractVGT, ExtractVITO or 
#ExtractGIMMS if possible

#local variables
ndvidirectory=paste(system.file("extdata/VITO_Mzimba",
    package="ndvits"), "/", sep="")
region="Mzimba"
Ystart=2004
Yend=2009
list="MzimbaList.txt"
period = 36
shape="SLP_Mzimba"
shapedir=paste(system.file("extdata/shape", package="ndvits"),
    "/", sep="")

#executing tolist to get a list
tolist(ndvidirectory, region, Ystart, Yend, 
    , outfile=list) 

#executing ExtractFiles
\donttest{signal = ExtractFile(shape, shapedir, list, outfile="TS.txt",
    period = period)}
}
\keyword{ts extraction}
