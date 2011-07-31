\name{sdNA}
\alias{sdNA}
\title{
sd function discarding NA values
}
\description{
Compute the standard deviation of a vector, discarding NA values.
}
\usage{
sdNA(x)
}
\arguments{
  \item{x}{vector}
}
\value{Return the standard deviation of \verb{x}}
\author{Romain Frelat, Bruno Gerard}
\seealso{
\code{\link{minNA}},\code{\link{meanNA}},\code{\link{maxNA}},\code{\link{minNA2}},\code{\link{maxNA2}}
}
\examples{
a=c(15,35,NA)
sdNA(a)
}
\keyword{tools}
