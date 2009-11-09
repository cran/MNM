\name{mv.ind.test}
\alias{mv.ind.test}
\Rdversion{1.1}

\title{
Independence Test
}
\description{
Tests for independence of two vectors using different scores.
}
\usage{
mv.ind.test(X, Y, score = "identity", method = "approximation", 
            n.simu = 1000, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric data frame or matrix. Must have the same number of rows as \code{Y}.
}
  \item{Y}{
a numeric data frame or matrix. Must have the same number of rows as \code{X}.
}
  \item{score}{
the score to be used. Possible are \code{identity}, \code{sign}, \code{symm} and \code{rank}.
}
  \item{method}{
method for for computation of the p-value for the spatial sign and spatial signed-rank tests.
                Possible are \code{approximation} and \code{permutation}.
}
  \item{n.simu}{
number of permutations if  \code{method="permutation"}.
}
  \item{na.action}{
a function which indicates what should happen when the data
          contain 'NA's. Default is to fail. 
}
}
\details{
This implements the independence tests as described in chapter 10 of the MNM book. Note that only inner test versions are implemented and that for the symmetrized sign score only
the approximative method for the computation of the p-value is available. 
}
\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the test statistic.}
  \item{parameter}{the degrees of freedom for the test statistic or the number of replications in the simulation.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified null hypothesis value.}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the two data matrices.}
}

\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } }

\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }

\examples{
X <- rmvt(150,diag(1,3),df=3)
Y <- rmvt(150, matrix(c(1,0.5,0.5,1),nrow=2),df=3)

mv.ind.test(X, Y)
mv.ind.test(X, Y, method = "p")

mv.ind.test(X, Y, score = "si")
mv.ind.test(X, Y, score = "si", method = "p")

mv.ind.test(X, Y, score = "r")
mv.ind.test(X, Y, score = "r", method = "p")

mv.ind.test(X, Y, score = "sy")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
