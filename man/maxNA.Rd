\name{maxNA}
\alias{maxNA}
\title{
max function discarding NA values
}
\description{
Compute the maximum of a vector, discarding NA values.
}
\usage{
maxNA(x)
}
\arguments{
  \item{x}{vector}
}
\value{Return the maximum of \verb{x}}
\author{Romain Frelat, Bruno Gerard}
\seealso{
\code{\link{maxNA2}},\code{\link{minNA}},\code{\link{meanNA}},\code{\link{sdNA}},\code{\link{minNA2}}
}
\examples{
a=c(15,35,NA)
maxNA(a)
}
\keyword{tools}
