\name{maxNA2}
\alias{maxNA2}
\title{
max function discarding NA values(2)
}
\description{
Compute the maximum of a two dimension vector, discarding NA values.
}
\usage{
maxNA2(x)
}
\arguments{
  \item{x}{a two dimension vector}
}
\value{Return the maximum of \verb{x}}
\author{Romain Frelat, Bruno Gerard}
\seealso{
\code{\link{maxNA}},\code{\link{minNA}},\code{\link{meanNA}},\code{\link{sdNA}},\code{\link{minNA2}}
}
\examples{
a=c(15,NA)
maxNA(a)
}
\keyword{tools}
