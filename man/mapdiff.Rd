\name{mapdiff}
\alias{mapdiff}
\title{Compute the difference between two maps}
\description{Compute the difference between two maps.\cr
If the map of the standard deviation is given, the difference is divided by the standard deviation.
}
\usage{
mapdiff(mapin, mapref, mapsd = FALSE, outfile = FALSE)
}
\arguments{
  \item{mapin}{raster map}
  \item{mapref}{raster map 2 (the reference)}
  \item{mapsd}{raster map conaining the standard deviation values.\cr 
If \var{mapsd}=FALSE, the standard deviation is not taken into account while computing the difference. By default, \var{mapsd} = FALSE.}
  \item{outfile}{file where the difference map is saved. 
If \var{outfile}=FALSE, the map is not saved but only returned. By default, \var{outfile} = FALSE.}
}
\details{
If mapsd is defined :
\deqn{diff = \frac{mapin - mapref}{mapsd}}
else :
\deqn{diff = mapin - mapref}
}
\value{
Returns the difference map
}
\author{Romain Frelat and Bruno Gerard}
\examples{
ndvidirectory=paste(system.file("extdata/VITO_Mzimba",
    package="ndvits"), "/", sep="")
region="Mzimba"
filein=timetoMap(ndvidirectory,region, 2004,1,1)
filein2=timetoMap(ndvidirectory,region, 2005,1,1)

#load the maps
map04=readpartGDAL(filein)
map05=readpartGDAL(filein2)

#computing the difference of the 2 maps
diff=mapdiff(map04, map05)

#plot the difference map :
savemap(diff)
}
\keyword{spatial}
