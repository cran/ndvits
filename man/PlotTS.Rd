\name{PlotTS}
\alias{PlotTS}
\title{
The function plots the time series one under another.
}
\description{
Time series are plotted one under another with the same scale of time.
}
\usage{
PlotTS(TS, outfile = FALSE, Ystart, period, title = "NDVI time series",
    nb = NULL)
}
\arguments{
  \item{TS}{Time serie or matrix of time serie containing only NDVI data (can be a single time series or multiple). }
  \item{outfile}{file where to save the graph. If FALSE (value by default), plot are shown in R.}
  \item{Ystart}{start year of the time series.}
  \item{period}{number of observations per year (for vito: period=36, for gimms: period=24).}
  \item{title}{title of the graphs. By default, title="NDVI time series".}
  \item{nb}{number of time series to plot in the same graph.\cr
By default, \var{nb} = NULL, it will be automatically chosen.}
}
\author{
Romain Frelat and Bruno Gerard
}
\seealso{
\code{\link{TimeSeriesAnalysis}},\code{\link{OverlayTS}}
}
\examples{
#loading the data
data(SLPSAs_ts)

#PlotTS
PlotTS(SLPSAs_ts$ts, Ystart=SLPSAs_ts$Ystart, 
    period=SLPSAs_ts$period, nb=6)
}
\keyword{plot}
