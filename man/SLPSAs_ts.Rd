\name{SLPSAs_ts}
\alias{SLPSAs_ts}
\docType{data}
\title{
NDVI time series extracted from Spot Vegetation data.
}
\description{
NDVI time series extracted from Spot Vegetation data.
The sites of interest are the South Asian study sites of CGIAR
Systemwide Livestock Programme. \cr
A kml file containing the location of the sites can be found at the following path : \cr
\code{print(system.file("extdata/shape/SLP-SAs.kml",package="ndvits"))}
}
\usage{data(SLPSAs_ts)}
\format{
  SLPSAs_ts is a list containing the 3 following objects :
\itemize{
\item ts : data frame containing the ndvi time series, one row per time series.\cr
\verb{rownames(ts)} contains the name of the site and \verb{colnames(ts)} contains the time.
\item period indicates the number of observations per year (period = 36).
\item Ystart indicates the starting year of the analysis (Ystart = 1999).
}
}
\source{
Spot Vegetation data comes from :
http://free.vgt.vito.be/
}
\examples{
data(SLPSAs_ts)

PlotTS(SLPSAs_ts$ts, Ystart=SLPSAs_ts$Ystart, 
    period=SLPSAs_ts$period, nb=6)

OverlayTS(SLPSAs_ts$ts, Ystart=SLPSAs_ts$Ystart, 
    period=SLPSAs_ts$period)

outfile = "SLPSAs-Metrics.txt"
outgraph = "SLPSAs-Metrics.pdf"
metrics = AnPhenoMetrics(SLPSAs_ts$ts, outfile, outgraph, 
    Ystart=SLPSAs_ts$Ystart, period=SLPSAs_ts$period)

#Be carefull, BFAST Decomposition takes a long time
#BfastDec(SLPSAs_ts$ts, Ystart=SLPSAs_ts$Ystart, 
#    period=SLPSAs_ts$period)
}
\keyword{datasets}
