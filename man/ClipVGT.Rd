\name{ClipVGT}
\alias{ClipVGT}
\title{
Clip Vegetation data
}
\description{
Clip the original map over a define rectangle.
This operation reduce the size of the map and the time to load it.
}
\usage{
ClipVGT(ndvidirectory, region, xlim, ylim, Ystart, Yend, outdirectory, regionout)
}
\arguments{
  \item{ndvidirectory}{directory of the original ndvi files}
  \item{region}{name of the region (as coded in the name of the file}
  \item{xlim}{longitude of the desired area}
  \item{ylim}{latitude of the desired area}
  \item{Ystart}{year of start of the clipping process}
  \item{Yend}{year of start of the clipping process}
  \item{outdirectory}{directory where to save the clipped maps}
  \item{regionout}{name of the region of the clipped area}
}
\details{
This function takes, as input, maps from FreeVGT (type=VITO_VGT) and save maps with "our" nomenclature (type=VITO_CLIP).\cr

FreeVGT VITO's nomenclature : (type=VITO_VGT)\cr
\verb{NDV_[YYYYMMDD]_[REGION]_Extract.tif}
\itemize{
\item{YYYY}{corresponging to the year (four digits)}
\item{MM}{corresponging to the months (two digits)}
\item{DD}{corresponging to the date of the period : "01", "11" or "21" (10-day mesures).}
\item{REGION}{the region selected while extracting maps}
}
Our nomenclature : : (type=VITO_CLIP)\cr
After a process of clipping the images in a particular region, we rename the files as following :\cr
\verb{[REGION][YYYY]M[MM]P[D].tif}
\itemize{
\item{REGION}{the region as entered while clipping data.}
\item{YYYY}{corresponging to the year (four digits)}
\item{MM}{corresponging to the months (two digits)}
\item{D}{corresponding to the period (0, 1 or 2).}
}
}
\value{Nothing is return, clipped map are saved in the \var{outdirectory}.}
\author{Romain Frelat and Bruno Gerard}
\seealso{
\code{\link{ExtractVGT}}, \code{\link{ExtractVito}}
}
\keyword{tools}
