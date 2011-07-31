\name{SLPSAs_full}
\alias{SLPSAs_full}
\docType{data}
\title{
Gross NDVI time series extracted from Spot Vegetation data.
}
\description{
Gross NDVI time series extracted from Spot Vegetation data.
The sites of interest are the South Asian study sites of CGIAR
Systemwide Livestock Programme. \cr
A kml file containing the location of the sites can be found at the following path : \cr
\code{print(system.file("extdata/shape/SLP-SAs.kml",package="ndvits"))}
}
\usage{data(SLPSAs_full)}
\format{
  SLPSAs_full is a list containing the 4 following objects :
\itemize{
\item full : a data.frame with 79 obs. of  397 variables \cr
\code{name} is a factor with levels \code{Assandh, India}, \code{Gharaunda}, \code{Goganda,India}, \code{Khansama, Bangladesh}, \code{Phulbari, Bangladesh}, \code{Sarada,India}.\cr
All the other fields are the time series as extracted directly from the NDVI images.
\item period indicates the number of observations per year (period = 36).
\item Ystart indicates the starting year of the analysis (Ystart = 1999).
\item maxNDVI indicates way of coding the ndvi images (maxNDVI = 255).
}
}
\source{
Spot Vegetation data comes from :
http://free.vgt.vito.be/
}
\examples{
data(SLPSAs_full)

ts = normNDVI(SLPSAs_full$full, SLPSAs_full$maxNDVI)

signal=STLperArea(ts, area = SLPSAs_full$full$name,
    Ystart = SLPSAs_full$Ystart, period = SLPSAs_full$period)
}
\keyword{datasets}
