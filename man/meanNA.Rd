\name{meanNA}
\alias{meanNA}
\title{
Mean function discarding NA values
}
\description{
Compute the mean of a vector, discarding NA values.
}
\usage{
meanNA(x)
}
\arguments{
  \item{x}{numerical vector}
}
\value{Return the mean of \verb{x}}
\author{Romain Frelat, Bruno Gerard}
\seealso{
\code{\link{minNA}},\code{\link{maxNA}},\code{\link{sdNA}},\code{\link{minNA2}},\code{\link{maxNA2}}
}
\examples{
a=c(15,35,NA)
meanNA(a)
}
\keyword{tools}
