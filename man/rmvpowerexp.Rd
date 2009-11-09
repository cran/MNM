\name{rmvpowerexp}
\alias{rmvpowerexp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Samples From a Power Exponential Distributions }
\description{
  Function to obtain random samples from a multivariate power exponential distribution. 
}
\usage{
rmvpowerexp(n, Location = rep(0, nrow(Scatter)), 
            Scatter = diag(length(Location)), Beta = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of random samples. }
  \item{Location}{Location vector of the  distribution.}
  \item{Scatter}{Scatter matrix of the distribution.}
  \item{Beta}{ shape parameter of the distribution.}
}
\details{
  The power exponential distribution is an elliptical distribution which can have light or heavy tails.
  \code{Beta = 1} yields a multivariate normal distribution, \code{Beta = 0.5} the multivariate Laplace distribution and
  with increasing \code{Beta} converges to a multivariate uniform distribution.
}
\value{
  a matrix.
}
\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } }
\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }

\seealso{ \code{\link[mvtnorm]{rmvnorm}}, \code{\link[mvtnorm]{rmvt}} }
\examples{
X1 <- rmvpowerexp(100,c(0,0,0),Beta = 0.5)
pairs(X1)
X2 <- rmvpowerexp(100,c(0,0,0),Beta = 1)
pairs(X2)
X3 <- rmvpowerexp(100,c(0,0,0),Beta = 10)
pairs(X3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ distribution }
\keyword{ multivariate}

