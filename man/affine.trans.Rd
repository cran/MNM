\name{affine.trans}
\alias{affine.trans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function For Affine Data Transformation }
\description{
  Function for transformations of the form \eqn{Ax + b}{Ax+b} or \eqn{A^{1/2}x + b}{A^(1/2)x+b.}
}
\usage{
affine.trans(X, A = diag(1, dim(X)[2]), b = rep(0, dim(X)[2]), 
             A.sqrt = FALSE, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data frame or matrix with p columns. }
  \item{A}{full rank p times p matrix. }
  \item{b}{numeric vector of length p. }
  \item{A.sqrt}{logical. If \code{TRUE} the symmetric square root of A will be used.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail. }
}

\value{
  a matrix.
}
\author{Klaus Nordhausen}
\examples{
data(iris)
IRIS <- iris[,1:4]
colMeans(IRIS)
cov(IRIS)
IRIS.trans <- affine.trans(IRIS, solve(cov(IRIS)), colMeans(IRIS),TRUE)
colMeans(IRIS.trans)
cov(IRIS.trans)
}
\keyword{ multivariate }

