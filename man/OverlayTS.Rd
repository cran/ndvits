\name{OverlayTS}
\alias{OverlayTS}
\title{Overlay time series}
\description{
Make one graph with all time series overlayed.
}
\usage{
OverlayTS(TS, Ystart, period = 36, title = "NDVI Time Series", outfile = FALSE)
}
\arguments{
  \item{TS}{matrix of time serie containing multiple NDVI time series.}
  \item{Ystart}{start year of the time series.}
  \item{period}{number of observations per year (for vito: period=36, for gimms: period=24).}
  \item{title}{title of the graphs. By default, title="NDVI time series".}
  \item{outfile}{file where to save the graph. If FALSE (value by default), plot are shown in R.}
}
\value{
This function makes a graph and returns nothing.
}
\author{
Romain Frelat and Bruno Gerard
}
\seealso{
\code{\link{PlotTS}}
}
\examples{
#loading the data
data(SLPSAs_ts)

#OverlayTS
OverlayTS(SLPSAs_ts$ts, Ystart=SLPSAs_ts$Ystart, 
    period=SLPSAs_ts$period)
}
\keyword{plot}
