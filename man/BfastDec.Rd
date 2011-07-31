\name{BfastDec}
\alias{BfastDec}
\title{BFAST decomposition}
\description{
Breaks for Additive Seasonal and Trend
}
\usage{
BfastDec(TS, Ystart, period, param, outfile = FALSE, ...)
}
\arguments{
  \item{TS}{matrix or data.frame containing the time series (one row per time series).
}
  \item{Ystart}{
starting year of the analysis (in four-digit format).
}
  \item{period}{
number of observations per year. \cr
For VGT's dataset : \var{period}=36, GIMMS's dataset: \var{period}=24.
}
  \item{param}{
list of parameters for Bfast decomposition
$season : "harmonic" or "dummy"
By default, param$season = "harmonic".
$h 
By default, param$h = period * 2/(length(as.numeric(TS[1, ]))).
$maxi 
By default, param$maxi = 10.
}
  \item{outfile}{
name of the file where the decomposition graphs are saved (extension ".pdf").\cr
If \var{outfile} = FALSE, the graphs are shown in R. By default, \var{outfile}=FALSE.
}
  \item{\dots}{further arguments sent to \verb{Bfast} function.}
}
\value{
res$seas
res$trend
}
\references{
Verbesselt, J., Hyndman, R., Newnham, G. and Culvenor, D. 2010 \emph{Remote Sensing of Environment}, Detecting trend and seasonal changes in satellite image time series, \bold{114}, 106--115

Verbesselt, J., Hyndman, R., Zeileis, A. and Culvenor, D. 2010 \emph{Remote Sensing of Environment}, Phenological change detection while accounting for abrupt and gradual trends in satellite image time series, \bold{114}, 2970--2980
}
\author{
Romain Frelat and Bruno Gerard
}
\seealso{
\code{\link{STLperArea}}
}
\examples{
data(SLPSAs_ts)

#Be carefull, BFAST Decomposition takes a long time
#BfastDec(SLPSAs_ts$ts, Ystart=SLPSAs_ts$Ystart, 
#    period=SLPSAs_ts$period)
}
\keyword{ts decomposition}
