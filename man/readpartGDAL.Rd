\name{readpartGDAL}
\alias{readpartGDAL}
\title{partial readGDAL}
\description{
Load geotiff maps partially, only the area of interest delimited by \var{xlim} and \var{ylim}.}
\usage{readpartGDAL(x, xlim = NULL, ylim = NULL, ...)}
\arguments{
  \item{x}{full path to the map to be loaded}
  \item{xlim}{
a two dimensions vector containing the longitude limits of the area of interest.\cr
\var{xlim}=NULL indicates no longitude limits. By default, \var{xlim}=NULL.\cr
Please note that \var{xlim} have to be in the same units and projection as the ndvi maps.
}
  \item{ylim}{
a two dimensions vector containing the latitude limits of the area of interest.\cr
\var{ylim}=NULL indicates no longitude limits. By default, \var{ylim}=NULL.\cr
Please note that \var{ylim} have to be in the same units and projection as the ndvi maps.
}
  \item{\dots}{further arguments sent to \code{\link{readGDAL}}.
}
}
\value{
Returns the map of the region of interest.
}
\author{
Romain Frelat and Bruno Gerard
}
\examples{
library(ndvits)
# Local variables
ndvidirectory=paste(system.file("extdata/VITO_Mzimba", package="ndvits"), "/", sep="")
region="Mzimba"
year=2004
#Define the limit of the area of interest :
xlim=c(33.5, 33.6)
ylim=c(-11.6, -11.4)
#call timetoMap to get the full path to the image
filein=timetoMap(ndvidirectory,region, year, 1, 1)

map = readpartGDAL(filein, xlim, ylim)

#display the loaded map
savemap(map)
}
\keyword{spatial}
