\name{STLperArea}
\alias{STLperArea}
\title{Seasonal Decomposition per Area}
\description{
Make a Seasonal Decomposition of Time Series by Loess per group of pixels
}
\usage{
STLperArea(TS, area, outfile = FALSE, Ystart, period = 36, fct = "mean", SGfilter = TRUE, nSG = "5,5", DSG = 0)
}
\arguments{
  \item{TS}{matrix or data.frame containing the time series (one row per time series).}
  \item{area}{factor which group pixels of the same area together (applying \verb{fct} per \verb{area})}
  \item{outfile}{name of the file where the decomposition is saved (extension ".pdf").\cr
If \var{outfile}=FALSE, the graphs are not saved but displayed in the R graphic console. By default, \var{outfile} = FALSE.}
  \item{Ystart}{starting year of the analysis (in four-digit format).}
  \item{period}{number of observations per year. \cr
For VGT's dataset : \var{period}=36, GIMMS's dataset: \var{period}=24.}
  \item{fct}{function to be applied to group points from the same \verb{area} : "mean", "max" or "min". By default \var{fct}="mean".}
  \item{SGfilter}{logical value indicating the application of Savitzky-Golay (SG) filter on time series
By default SGfilter=TRUE}
  \item{nSG}{parameters of the windows of the SG filter.\cr
By default nSG="5,5" (windows 5,5).}
  \item{DSG}{parameters of the degree polynome of the SG filter. 
By default, DSG=0.}
}
\value{
The output of STLperArea is the time series computed (mean/max of the time series grouped by \verb{area} and with/without Savitzky-Golay filter).
}
\author{
Romain Frelat and Bruno Gerard
}
\seealso{
\code{\link{TimeSeriesAnalysis}},\code{\link{BfastDec}}
}
\examples{
#loading the data
data(SLPSAs_full)

#formating the time series
ts = normNDVI(SLPSAs_full$full, SLPSAs_full$maxNDVI)

#STLperArea
signal=STLperArea(ts, area = SLPSAs_full$full$name, 
    Ystart = SLPSAs_full$Ystart, period = SLPSAs_full$period)
}
\keyword{ts decomposition}
