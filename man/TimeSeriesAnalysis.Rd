\name{TimeSeriesAnalysis}
\alias{TimeSeriesAnalysis}
\title{NDVI Time Series Extraction and Analysis}
\description{
Extract NDVI time series from different data source.
Group points or polygon together, compute savitzky-Golay filter.
Also exectute a basic analysis (STL decomposition).
}
\usage{
TimeSeriesAnalysis(shapefile, shapedir, ndvidirectory, region,
    Ystart, Yend, outfile = "TS.txt", outfile2 = "TS.pdf", 
    outfile3 = FALSE, shapeext = "shp", fct = "mean", 
    SGfilter = TRUE, nSG = "5,5", DSG = 0, 
    title = "NDVI time series", type = "VITO_CLIP", nb = NULL)
}
\arguments{
  \item{shapefile}{
name of the shape/kml file with multi-points or polygons from where time series are extracted. \cr
For a shapefile called \file{name.shp} , \var{shapefile}="name".\cr
For a kml file, write the name of the layer.
}
  \item{shapedir}{
path to the shape/kml file. \cr
For a shapefile, \var{shapedir} is the path to the folder containing the shapefile. 
For example if the directory of the shapefile is \file{C:/Dir/name.shp}, \var{shapefile}="C:/Dir/". \cr
For a kml file, \var{shapedir} is the full directory of the kml file. For example \var{shapefile}="C:/Dir/name.kml". \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and "/" in UNIX OS.}
  \item{ndvidirectory}{
path to the ndvi images folder (GeoTIFF files). \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and with "/" in UNIX OS
}
  \item{region}{
name of the region, same as in the name of the ndvi images. \cr
Please refer to details section for more information on the nomenclature.
}
  \item{Ystart}{
starting year of the analysis (in four-digit format).
}
  \item{Yend}{
ending year of the analysis (in four-digit format).
}
  \item{outfile}{
name of the file where the time series are saved (extension ".txt").\cr By default, \var{outfile} = " .".
}
  \item{outfile2}{
file where the seasonal decomposotion of the time serie is saved (extenstion ".pdf"). By default, outfile2="TS.pdf". }
  \item{outfile3}{
 file where the graphs of the different time serie (stocked) is saved (extenstion ".pdf").
If outfile3=FALSE, graphs are displayed in R (but not saved). By default, outfile3=FALSE
}
  \item{shapeext}{
extension of the shape/kml file : "shp" for a Esri shapefile or "kml" for a kml file from Google earth .\cr 
By default, \var{shapeext} = "shp".
}
  \item{fct}{
function applied while grouping points of a same categorie together : "mean","max" or "min".
By default, fct="mean".
}
  \item{SGfilter}{
TRUE to apply a Savitzky-Golay filter to the time series, FALSE otherwize.
By default, SGfilter=TRUE.
}
  \item{nSG}{
moving window of the Savitzky-Golay filter.
By default, nSG="5,5" (symetric windows taking into account 5 points backward and 5 points forward.
}
  \item{DSG}{
polynomial degree of the function use in the Savitzky-Golay filter.
By default, DSG=0.
}
  \item{title}{
title of the graph (saved in \var{outfile3})
By default, title="NDVI time series".
}
  \item{type}{
type of nomenclature followed by the NDVI images : "VITO_CLIP", "VITO_VGT" or "GIMMS" (see details for more information).\cr
By default, \var{type} = "VITO_CLIP".
type of nomenclature followed by the NDVI images : "GIMMS", "VITO_VGT", "VITO_CLIP", "TEXT" or "FILES" (see details for more information).\cr
By default, \var{type} = "VITO_CLIP".
}
  \item{nb}{
number of time series plotted in the same graph by the function \code{\link{PlotTS}} saved in the file \var{outfile3}. 
By default, \var{nb} = NULL : will be automatically chosen depending on the number of time series.
}
}
\details{
More about the variable \var{type} : \cr
\var{type}=GIMMS : nomenclature from Global Land Cover Facility \cr
\verb{[RE][YY][MMM][DDD].[SAT]-VIg.tif }
\itemize{
\item{RE}{identifies the continent (AF for Africa,AZ for Australia and New Zealand, EA for Eurasia, NA for North America and SA for South America and Central America)}
\item{YY}{the two last digits of the year}
\item{MMM}{three-letter abbreviation of the month}
\item{DDD}{identifies the composite period, 15a denotes the days 1-15 of the month, and 15b denotes the days from 16 to the end of the month.}
\item{SAT}{identifies the satellite, NOAA-7, 9, 11, 14, 16, or 17 from which the data originated.\cr}
}
\var{type}=VITO_VGT : FreeVGT VITO's nomenclature \cr
\verb{NDV_[YYYYMMDD]_[REGION]_Extract.tif }
\itemize{
\item{YYYY}{identifies the year (four digits)}
\item{MM}{identifies the months (two digits)}
\item{DD}{identifies the date of the period : "01", "11" or "21" (10-day composite).}
\item{REGION}{idendifies the region selected while extracting maps\cr}
}
\var{type}=VITO_CLIP : nomenclature after clipping the images in a particular region.\cr
\verb{[REGION][YYYY]M[MM]P[D].tif }
\itemize{
\item{REGION}{the region as entered while clipping data.}
\item{YYYY}{identifies the year (four digits)}
\item{MM}{identifies the months (two digits)}
\item{D}{identifies the period (0, 1 or 2).\cr}
}
\var{type} = FILES : 
\var{ndvidirectory} should contain a file respecting the following format :\cr
\verb{N}\cr
\verb{path\image1.tif}\cr
\verb{path\image2.tif}\cr
\verb{...}\cr
\verb{path\imageN.tif}\cr
with \verb{N} the number of image files. \cr
The list of images should be in chronological order, and include the complete path to the file, 
unless the files are located in the working directory from which R is run. Only the N first images in the file will be processed.
nomenclature after clipping the images in a particular region.\cr
}
\value{
data.frame with the extracted and processed time series.
}
\author{
R. Frelat, B. Gerard
}
\seealso{
\code{\link{ExtractFile}}, \code{\link{ExtractGIMMS}}, \code{\link{ExtractVGT}}, \code{\link{ExtractVito}}, 
\code{\link{normNDVI}}, \code{\link{PlotTS}}, \code{\link{STLperArea}}
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
}
\keyword{ts extraction}
