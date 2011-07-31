\name{ExtractGIMMS}
\alias{ExtractGIMMS}
\title{Extract time series from GIMMS data}
\description{Extract time series from GIMMS data.}
\usage{ExtractGIMMS(shapefile, shapedir, ndvidirectory, region,
    outfile, Ystart, Yend, shapeext = "shp")}
\arguments{
  \item{shapefile}{name of the shape/kml file with multi-points or polygons from where time series are extracted. \cr
For a shapefile called \file{name.shp} , \var{shapefile}="name".\cr
For a kml file, write the name of the layer.}
  \item{shapedir}{path to the shape/kml file. \cr
For a shapefile, \var{shapedir} is the path to the folder containing the shapefile. 
For example if the directory of the shapefile is \file{C:/Dir/name.shp}, \var{shapefile}="C:/Dir/". \cr
For a kml file, \var{shapedir} is the full directory of the kml file. For example \var{shapefile}="C:/Dir/name.kml". \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and "/" in UNIX OS.
}
  \item{ndvidirectory}{
path to the GIMMS ndvi images folder (GeoTIFF files). \cr
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
GIMMS nomenclature (from Global Land Cover Facility) : (type=GIMMS)\cr
\verb{[RE][YY][MMM][DDD].[SAT]-VIg.tif}
\itemize{
\item{RE}{identifies the continent (AF for Africa,AZ for Australia and New Zealand, EA for Eurasia, NA for North America and SA for South America and Central America)}
\item{YY}{corresponging to the year (two last digits)}
\item{MMM}{abbreviation of the month (three letters)}
\item{DDD}{identifies the composite period, 15a denotes the days 1-15 of the month, and 15b denotes the days from 16 to the end of the month.}
\item{SAT}{identifies the satellite, NOAA-7, 9, 11, 14, 16, or 17 from which the data originated.}
}
}
\value{Return the time series of the points/polygons in the shapefile over the desired period.}
\author{Romain Frelat and Bruno Gerard}
\seealso{
\code{\link{TimeSeriesAnalysis}}
}
\keyword{ts extraction}
