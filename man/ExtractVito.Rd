\name{ExtractVito}
\alias{ExtractVito}
\title{
Extract time series from clipped VITO data
}
\description{
Extract time series from clipped VITO data.
}
\usage{
ExtractVito(shapefile, shapedir, ndvidirectory, region, 
    outfile = "TS.txt", Ystart, Yend, shapeext = "shp")
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
  \item{ndvidirectory}{path to the ndvi images folder (GeoTIFF files). \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and with "/" in UNIX OS.}
  \item{region}{name of the region, same as in the name of the ndvi images. \cr
Please refer to details section for more information on the nomenclature.}
  \item{outfile}{name of the file where the time series are saved (extension ".txt").\cr By default, \var{outfile} = " .".}
  \item{Ystart}{starting year of the analysis (in four-digit format).}
  \item{Yend}{ending year of the analysis (in four-digit format).}
  \item{shapeext}{extension of the shape/kml file : "shp" for a Esri shapefile or "kml" for a kml file from Google earth .\cr 
By default, \var{shapeext} = "shp".}
}
\details{Files should follow this nomenclature (type=VITO_CLIP):\cr
After a process of clipping the images in a particular region, we rename the files as following :\cr
\verb{[REGION][YYYY]M[MM]P[D].tif}
\itemize{
\item{REGION}{the region as entered while clipping data.}
\item{YYYY}{corresponging to the year (four digits)}
\item{MM}{corresponging to the months (two digits)}
\item{D}{corresponding to the period (0, 1 or 2).}
}
}
\value{Return the time series of the points/polygons in the shapefile over the desired period.}
\author{Romain Frelat and Bruno Gerard}
\seealso{
\code{\link{TimeSeriesAnalysis}}, \code{\link{ClipVGT}}
}
\examples{
#local variables
ndvidirectory=paste(system.file("extdata/VITO_Mzimba",
    package="ndvits"), "/", sep="")
region="Mzimba"
Ystart=2004
Yend=2009
shape="SLP_Mzimba"
shapedir=paste(system.file("extdata/shape", package="ndvits"),
    "/", sep="")

#executing ExtractFiles
signal = ExtractVito(shape, shapedir, ndvidirectory, region, 
    outfile = "TS-Mzimba.txt", Ystart = Ystart, Yend = Yend)
}
\keyword{ts extraction}
