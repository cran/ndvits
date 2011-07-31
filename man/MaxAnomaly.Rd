\name{MaxAnomaly}
\alias{MaxAnomaly}
\title{Maximum Anomaly}
\description{
The function MaxAnomaly compute the anomaly map of the difference between the maximum of each year and their mean.
}
\usage{
MaxAnomaly(ndvidirectory, region, Ystart, Yend, 
    outname = "anomaly", outext = "show", xlim = NULL, 
    ylim = NULL, type="VITO_CLIP", shapefile = NULL, 
    shapedir = NULL, shapeext = "shp", label = FALSE, 
    pal = "Spectral")
}
\arguments{
  \item{ndvidirectory}{path to the ndvi images folder (GeoTIFF files). \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and with "/" in UNIX OS.}
  \item{region}{name of the region, same as in the name of the ndvi images. \cr
Please refer to details section for more information on the nomenclature.}
  \item{Ystart}{starting year of the analysis (in four-digit format).}
  \item{Yend}{ending year of the analysis (in four-digit format).}
  \item{outname}{prefix of the name of the files where the anomaly maps are saved.\cr
By default \var{outname}="anomaly"}
  \item{outext}{extension of the created maps. \cr
If \var{outext} = "show", maps are displayed in the R graphics console.\cr
Choose between "jpg", "png" or "pdf" to save the map in the corresponding format.\cr
By default, \var{outext} = "show".}
  \item{xlim}{a two dimensions vector containing the longitude limits of the area of interest.\cr
\var{xlim}=NULL indicates no longitude limits. By default, \var{xlim}=NULL.\cr
Please note that \var{xlim} have to be in the same units and projection as the ndvi maps.}
  \item{ylim}{a two dimensions vector containing the latitude limits of the area of interest.\cr
\var{ylim}=NULL indicates no longitude limits. By default, \var{ylim}=NULL.\cr
Please note that \var{ylim} have to be in the same units and projection as the ndvi maps.}
  \item{type}{type of nomenclature followed by the NDVI images : "VITO_CLIP", "VITO_VGT" or "GIMMS" (see details for more information).\cr
By default, \var{type} = "VITO_CLIP".}
  \item{shapefile}{name of the shape/kml file with objects (points or polygons) to be displayed on the maps.\cr
For a shapefile called \file{name.shp} , \var{shapefile}="name".\cr
For a kml file, write the name of the layer.\cr
If \var{shapefile}=NULL, no object will be added to the map.\cr
By default, \var{shapefile}=NULL.}
  \item{shapedir}{path to the shape/kml file. \cr
For a shapefile, \var{shapedir} is the path to the folder containing the shapefile. 
For example if the directory of the shapefile is \file{C:/Dir/name.shp}, \var{shapefile}="C:/Dir/". \cr
For a kml file, \var{shapedir} is the full directory of the kml file. For example \var{shapefile}="C:/Dir/name.kml". \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and "/" in UNIX OS.}
  \item{shapeext}{extension of the shape/kml file : "shp" for a Esri shapefile or "kml" for a kml file from Google earth .\cr 
By default, \var{shapeext} = "shp".}
  \item{label}{name of the label of the points to be displayed on the map. \cr 
For a shapefile, please indicates the name of the category you want to see.
For a kml file, \var{label}="Name".
If label=FALSE, no label will be plotted. By default, \var{label}=FALSE.}
  \item{pal}{name of the palette of colors used to make the map from the package \pkg{RColorBrewer}.\cr
To see the different palettes available for your map, please type \code{display.brewer.all()}.\cr
By default, \var{pal}="Spectral".}
}
\details{
The map of the max of each year between \verb{Ystart} and \verb{Yend} is computed. Then, the mean is calculated.  The anomaly map comes from the difference between the mean of the max and the max of the specific year.  It is saved automatically in geotiff and displayed by default (show=TRUE) in the graphic of R.  To go from one graph to another, click on the graph or press enter in the R console.\cr

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
\author{Romain Frelat and Bruno Gerard}
\seealso{
\code{\link{maplocalstat}}, \code{\link{PeriodAnomaly}}
}
\examples{
#local variables
ndvidirectory=paste(system.file("extdata/VITO_Mzimba",
    package="ndvits"), "/", sep="")
region="Mzimba"
Ystart=2004
Yend=2009

#executing MaxAnomaly
\donttest{MaxAnomaly(ndvidirectory, region, Ystart, Yend,
    outname="Mzimba", outext="show")}

#adding points to the map
shape="SLP_Mzimba"
shapedir=paste(system.file("extdata/shape", package="ndvits"),
    "/", sep="")
\donttest{MaxAnomaly(ndvidirectory, region, Ystart, Yend,
    outname="Mzimba", outext="show", shapefile=shape,
    shapedir=shapedir, label="Village")}
}
\keyword{anomaly map}
