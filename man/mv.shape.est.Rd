\name{mv.shape.est}
\alias{mv.shape.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Shape Matrices }
\description{
  Shape matrix estimates using different score functions.
}
\usage{
mv.shape.est(X, score = "identity", estimate = "outer", 
             location = NULL, na.action = na.fail, ...)
}
\arguments{
  \item{X}{a numeric data frame or matrix. }
  \item{score}{score to be used. Can be either \code{"identity"}, \code{"sign"}, \code{"symmsign"} or \code{"rank"}.}
  \item{estimate}{can be \code{"outer"} or \code{"inner"}.}
  \item{location}{If NULL the location vector is estimated. Alternatively a numeric p vector of location. }
  \item{na.action}{ a function which indicates what should happen when
    the data contain 'NA's. Default is to fail. }
  \item{\dots}{further arguments passed to or from other methods. }
}
\details{
  This functions returns different shape matrices depending on the score function chosen. For details see chapter 9 of the MNM book.
}
\value{
 a matrix
}
\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } 

\cite{Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: The Package MNM, \emph{Journal of Statistical Software}, \bold{43}, 1-28.}
}
\author{Klaus Nordhausen}


\seealso{ \code{\link[stats]{cov}}, \code{\link[ICSNP]{tyler.shape}}, \code{\link[ICSNP]{duembgen.shape}}, \code{\link[ICSNP]{HR.Mest}}, \code{\link[SpatialNP]{spatial.shape}}}

\examples{
data(iris)
IRIS <- iris[,1:4]
mv.shape.est(IRIS, "sign")
mv.shape.est(IRIS, "symmsign", "o")
mv.shape.est(IRIS, "rank")
}
\keyword{multivariate}
\keyword{nonparametric}
