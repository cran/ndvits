\name{tolist}
\alias{tolist}
\title{To list}
\description{
Create a file containing a list with the directories of maps.
This file can be used in Timesat or with the function \code{\link{ExtractFile}}.
}
\usage{
tolist(ndvidirectory, region, Ystart, Yend, outfile = "list.txt", type = "VITO_CLIP")
}
\arguments{
  \item{ndvidirectory}{path to the ndvi images folder (GeoTIFF files). \cr
Please note that paths to files in R can be specified with either "/" or "\\\\" in Windows OS and with "/" in UNIX OS.}
  \item{region}{name of the region, same as in the name of the ndvi images. \cr
Please refer to details section for more information on the nomenclature.}
  \item{Ystart}{starting year of the analysis (in four-digit format).}
  \item{Yend}{ending year of the analysis (in four-digit format).}
  \item{outfile}{name of the file where the list is saved (extenstion ".txt"). By default, \var{outfile}="list.txt".}
  \item{type}{type of nomenclature followed by the NDVI images : "VITO_CLIP", "VITO_VGT" or "GIMMS" (see details for more information).\cr
By default, \var{type} = "VITO_CLIP".}
}
\details{
The format of the output file is as follows:\cr
\verb{N}\cr
\verb{path\image1.tif}\cr
\verb{path\image2.tif}\cr
\verb{...}\cr
\verb{path\imageN.tif}\cr
with \verb{N} the number of image files. \cr
File names of image are in chronological order.
It is important, for further utilisation, that \var{ndvidirectory} include the complete path to the file. \cr 

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
Nothing return.
}
\references{
Jonsson, P. and Eklundh, L. 2004, \emph{Computers & Geosciences}, TIMESAT a program for analyzing time-series of satellite sensor data, \bold{30}, 833--845
}
\author{Romain Frelat and Bruno Gerard}
\seealso{\code{\link{ExtractFile}}}
\examples{
#local variables
ndvidirectory=paste(system.file("extdata/VITO_Mzimba",
    package="ndvits"), "/", sep="")
region="Mzimba"
Ystart=2004
Yend=2009

#executing tolist
tolist(ndvidirectory, region, Ystart, Yend, 
    , outfile="MzimbaList.txt") 

}
\keyword{tools}
