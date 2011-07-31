\name{pointtobuffer}
\alias{pointtobuffer}
\title{
Point to buffer
}
\description{
Create a shapefile of points around the points given in the kml/shapefile.
}
\usage{
pointtobuffer(shapefile, shapedir, ndvidirectory, region, Ystart, shapeext="shp", outshape = "buffer", outdir = ".", rad = 1, type="VITO_CLIP")
}
\arguments{
  \item{shapefile}{name of the shape/kml file with multi-points from where the buffer is created. \cr
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
  \item{Ystart}{starting year of the analysis (in four-digit format).}
  \item{shapeext}{extension of the shape/kml file : "shp" for a Esri shapefile or "kml" for a kml file from Google earth .\cr 
By default, \var{shapeext} = "shp".}
  \item{outshape}{name of the new shapefile where the buffer is saved.\cr
By default, \var{outshape} = "buffer".}
  \item{outdir}{path to the folder where the new shapefile is saved.\cr
By default, \var{outdir} = ".".}
  \item{rad}{radius, maximum distance in pixels.\cr
New points are added to the shapefile if the distance to the points is smaller than \var{rad}.\cr
By default, \var{rad}=1, 8 points will be added around each points in the given shapefile.}
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
\item{D}{ : identifies the period (0, 1 or 2).\cr}
}
}
\value{
Nothing is return, a shapefile is created and saved in the folder \verb{outdir}, under the name \verb{outshape}.
}
\author{
Romain Frelat and Bruno Gerard
}
\seealso{
\code{\link{polytopoints}}
}
\examples{
#local variables
shape="SLP_Mzimba"
shapedir=paste(system.file("extdata/shape", package="ndvits"),
    "/", sep="")
ndvidirectory=paste(system.file("extdata/VITO_Mzimba",
    package="ndvits"), "/", sep="")
region="Mzimba"
Ystart=2004
newshape="BufferMzimba"
dirshape="."
#the shapefile will be saved in the current directory.

pointtobuffer(shape, shapedir, ndvidirectory, region, Ystart,
   outshape = newshape, outdir = dirshape) 
}
\keyword{spatial}
