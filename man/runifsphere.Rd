\name{runifsphere}
\alias{runifsphere}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Random Samples From the Unit Sphere}
\description{
  Function to sample uniformly distributed observations on the unit sphere.
}
\usage{
runifsphere(n, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{ number of random samples.}
  \item{p}{ dimension of the unit sphere. }
}

\value{
  a matrix.
}
\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } }
\author{Klaus Nordhausen}


\examples{
X <- runifsphere(100,2)
plot(X, pty = "s")
}
\keyword{ distribution }
\keyword{ multivariate}
