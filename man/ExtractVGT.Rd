\name{ExtractVGT}
\alias{ExtractVGT}
\title{
Extract time series from VITO Vegetation data
}
\description{
Extract time series from VITO Vegetation data
}
\usage{
ExtractVGT(shapefile, shapedir, ndvidirectory, region, 
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
  \item{ndvidirectory}{path to the SPOT VGT ndvi images folder (GeoTIFF files). \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and with "/" in UNIX OS.}
  \item{region}{name of the region, same as in the name of the ndvi images. \cr
Please refer to details section for more information on the nomenclature.}
  \item{outfile}{name of the file where the time series are saved (extension ".txt").\cr By default, \var{outfile} = " .".}
  \item{Ystart}{starting year of the analysis (in four-digit format).}
  \item{Yend}{ending year of the analysis (in four-digit format).}
  \item{shapeext}{extension of the shape/kml file : "shp" for a Esri shapefile or "kml" for a kml file from Google earth .\cr 
By default, \var{shapeext} = "shp".}
}
\details{
FreeVGT VITO's nomenclature : (type=VITO_VGT)\cr
\verb{NDV_[YYYYMMDD]_[REGION]_Extract.tif}
\itemize{
\item{YYYY}{corresponging to the year (four digits)}
\item{MM}{corresponging to the months (two digits)}
\item{DD}{corresponging to the date of the period : "01", "11" or "21" (10-day mesures).}
\item{REGION}{the region selected while extracting maps}
}
}
\value{
Return the time series of the points/polygons in the shapefile over the desired period.
}
\references{
VGT ref
}
\author{
Romain Frelat and Bruno Gerard
}
\seealso{
\code{\link{TimeSeriesAnalysis}}, \code{\link{ClipVGT}}
}
\keyword{ts extraction}
