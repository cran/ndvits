\name{periodtoMap}
\alias{periodtoMap}
\title{
period to full path of the NDVI map.
}
\description{
Returns the full path of the map of a given period and a given year.
}
\usage{
periodtoMap(ndvidirectory, region, year, period, type = "VITO_CLIP")
}
\arguments{
  \item{ndvidirectory}{path of the ndvi image folder (GeoTIFF files). /cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and with "/" in UNIX OS.}
  \item{region}{name of the region, same as in the name of the ndvi images.}
  \item{year}{year of the map to load}
  \item{period}{period of the map to load. Be careful, for GIMMS data, period is between 1 and 24 while for VITO data, period is between 1 and 36.}
  \item{type}{type of nomenclature followed by the NDVI images : "VITO_CLIP", "VITO_VGT" or "GIMMS" (see details for more information).\cr
By default, \var{type} = "VITO_CLIP".}
}
\details{
More about the variable \var{type} : \cr
\var{type} = GIMMS : nomenclature from Global Land Cover Facility \cr
\verb{[RE][YY][MMM][DDD].[SAT]-VIg.tif }
\itemize{
\item{RE}{ : identifies the continent (AF for Africa,AZ for Australia and New Zealand, EA for Eurasia, NA for North America and SA for South America and Central America)}
\item{YY}{ : the two last digits of the year}
\item{MMM}{ : three-letter abbreviation of the month}
\item{DDD}{ : identifies the composite period, 15a denotes the days 1-15 of the month, and 15b denotes the days from 16 to the end of the month.}
\item{SAT}{ : identifies the satellite, NOAA-7, 9, 11, 14, 16, or 17 from which the data originated.\cr}
}
\var{type} = VITO_VGT : FreeVGT VITO's nomenclature \cr
\verb{NDV_[YYYYMMDD]_[REGION]_Extract.tif }
\itemize{
\item{YYYY}{ : identifies the year (four digits)}
\item{MM}{ : identifies the months (two digits)}
\item{DD}{ : identifies the date of the period : "01", "11" or "21" (10-day composite).}
\item{REGION}{ : idendifies the region selected while extracting maps\cr}
}
\var{type} = VITO_CLIP : nomenclature after clipping the images in a particular region.\cr
\verb{[REGION][YYYY]M[MM]P[D].tif }
\itemize{
\item{REGION}{ : the region as entered while clipping data.}
\item{YYYY}{ : identifies the year (four digits)}
\item{MM}{ : identifies the months (two digits)}
\item{D}{ : identifies the period (0, 1 or 2).}
}
}
\value{
Returns the complete path of the desired map.
If the map doesn't exist, an error message is displayed.
}
\author{
Romain Frelat and Bruno Gerard
}
\seealso{
\code{\link{timetoMap}}
}
\examples{
library(ndvits)
# Local variables
ndvidirectory=paste(system.file("extdata/VITO_Mzimba", package="ndvits"), "/", sep="")
region="Mzimba"
year=2004
period=15

filein=periodtoMap(ndvidirectory,region, year, period)
print(filein)
}
\keyword{tools}
