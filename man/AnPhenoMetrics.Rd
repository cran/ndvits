\name{AnPhenoMetrics}
\alias{AnPhenoMetrics}
\title{
Extracting annual phenological metrics
}
\description{
This function extracts annual phenological metrics
}
\usage{
AnPhenoMetrics(TS, outfile, outgraph, Ystart, period, SOSth = 0.5, EOSth = 0.5)
}
\arguments{
  \item{TS}{matrix or data.frame containing the time series (one row per time series).}
  \item{outfile}{name of the file where the metrics are saved (extension ".txt"). }
  \item{outgraph}{name of the file where the graphs are saved (extension ".pdf").\cr}
  \item{Ystart}{starting year of the analysis (in four-digit format).}
  \item{period}{number of observations per year. \cr
For VGT's dataset : \var{period}=36, GIMMS's dataset: \var{period}=24.}
  \item{SOSth}{threshold for the detection of the start of the season. \cr By default SOSth=0.5}
  \item{EOSth}{threshold for the detection of the start of the season. \cr By default SOSth=0.5}
}
\details{
We proceed in 3 mean steps :
\enumerate{
 \item Detecting the number of minimums and maximums with the mean signal.\cr
The mean signal over the year is computed. Then we detect the minimums and the maximums.
A maximum is detected if it is a maximum over a window of 6 measures and if its value is higher than the mean of the mean signal.
Similarly, a minimum is detected if it is a minimum over a window of 6 measures and if its value is lower than the mean of the mean signal.
Then, a routine check if minimums and maximums are one after another and if there are the same number of maximums and minimums.

If the global mean of the time series is higher than 0.7 or smaller than 0.2, further calculations are stopped.
The signal is the one of the bare soil or forest and metrics became very difficult to compute, no more reliable.  

 \item Focusing on individual year time series to detect minimums and maximums.\cr
The main issue is to handle the full season whenever its starts and its stopped (very often up to the next year).
Our solution is to start the time series two months before the global minimum. 
If two minimums (two seasons), we take into account the minimum of them.
For each year, we consider a time series with a length of 16 months to be sure to capture the full season.
The yearly minimums/maximums are detected within a two months windows around the minimums/maximums of the mean signal (detected during the first step).

 \item Calculating the phenological metrics. \cr
Every year, four metrics are extracted from the time series :
\itemize{
 \item SOS : Start Of the growing Season \cr
We use a threshold method. The day of the start of the season is detected when the NDVI value has increase by $50\%$ of the distance between the "left" minimum and the maximum.
The value of the threshold can be customized by the user.
 \item EOS : End Of the growing Season \cr
Similarly than SOS, the day of the end of the season is detected when the NDVI value has decrease by $50\%$ of the difference between the "right" minimum and the maximum.
The value of the threshold can be customized by the user.
 \item LOS : Length of the growing Season, $LOS = EOS - SOS$.\cr
LOS a number of period.
 \item cumNDVI : cumulative NDVI.\cr
integral under the NDVI curve during the growing season.
There isn't a consensus on the way to compute the cumulative NDVI.
The minimum value has to be set.
For the moment, we compute the full integral under the curves.
Hopefully, soon there will be an option to consider either 0.1 as minimum (to delete bare soil effect) or 
the mean of the SOS and EOS NDVI value (to take into account only the curve of the growing season).
}
}
}
\value{
A text file is saved with phenological metrics \verb{outfile} and a list is returned with :
\itemize{
\item names : vector with the names of the sites (rownames of the following matrix).
\item year : vector with the year (colnames of the following matrix).
\item mld : date of the left minimum (number of the period).
\item M : maximum NDVI value over the considered period.
\item SOS : Start Of the Season (number of the period)
\item EOS : End Of the Season (number of the period)
\item LOS : Length Of the Season (number of the period)
\item cumNDVI : integral under the NDVI curve during the growing season.
}
}
\references{
Brown, M. E. and de Beurs K. and Vrieling A. 2010 \emph{Remote Sensing of Environment}, The response of African land surface phenology to large scale climate oscillations, \bold{114}, 2286--2296
}
\author{
Romain Frelat and Bruno Gerard
}

\examples{
#loading the data
data(SLPSAs_ts)

#defining local variables
outfile = "SLPSAs-Metrics.txt"
outgraph = "SLPSAs-Metrics.pdf"

#extracting metrics
metrics = AnPhenoMetrics(SLPSAs_ts$ts, outfile, outgraph, 
    Ystart=SLPSAs_ts$Ystart, period=SLPSAs_ts$period)
}
\keyword{metrics}
