\name{plot.mvloc}
\Rdversion{1.1}
\alias{plot.mvloc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plotting Method for mvloc Objects
}
\description{
Visually presents and compares different multivariate location estimates and their confidence ellipsoids.
}
\usage{
\method{plot}{mvloc}(x, est2 = NULL, est3 = NULL, X = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class \code{mvl1lm}.}
  \item{est2}{An optional additional location estimate. A list with the components \code{location}, \code{vcov} and \code{est.name}, for example an object of class 'mvloc'.}
  \item{est3}{An optional additional location estimate. A list with the components \code{location}, \code{vcov} and \code{est.name}, for example an object of class 'mvloc'.}
  \item{X}{a numeric data frame or matrix. Optional data points on which the estimates could have been based.}
  \item{\dots}{optional plotting arguments. For details see \code{\link{plotMvloc} }.}
}
\details{
  The figure can be used to compare different multivariate location estimates. The location of the legend is currently problematic and it is
  recommended that the user should provide the coordinates for the legend. The function calls \code{\link{plotMvloc}}. 
}
\author{Klaus Nordhausen}

\seealso{ \code{\link{plotMvloc} }}

\examples{
X <- rmvt(50, diag(c(1, 2)), 3)

est1 <- mv.1sample.est(X)
est2 <- mv.1sample.est(X, "sign")
est3 <- mv.1sample.est(X, "rank", "inner")

plot(est1)
plot(est1, est2, est3, X, alim="b", lty.ell=1:3, pch.ell=14:16)
plot(est1, est2, est3, X, alim="e")
}

\keyword{methods}
