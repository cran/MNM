\name{mv.Csample.test}
\alias{mv.Csample.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{C Sample Test of Location }
\description{
  Several samples location tests using different scores.
}
\usage{
mv.Csample.test(X, g, score = "identity", stand = "outer", 
                method = "approximation", n.simu = 1000, 
                na.action = na.fail, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data frame or matrix of response values. }
  \item{g}{a factor with at least two levels.}
  \item{score}{the score to be used. Possible choices are \code{identity}, \code{sign} and \code{rank}.}
  \item{stand}{ the standardization method used. Possible choices are \code{outer} and \code{inner}.}
  \item{method}{method for the computation of the p-value for the spatial sign and spatial signed-rank tests.
                Possible choices are \code{approximation} and \code{permutation}.}
  \item{n.simu}{number of simulated sign changes if  \code{method="permutation"}.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail. } 
  \item{\dots}{arguments that can be passed on to functions used for the estimation of the spatial signs and spatial ranks. }
}
\details{
  This implements the location tests based on identity, sign or rank scores as described in chapter 11 of the MNM book.
}
\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the test statistic.}
  \item{parameter}{the degrees of freedom for the test statistic or the number of replications in the simulation.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified hypothesized value of the location.}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data set and of the grouping vector.}
}
\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } 

\cite{Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: The Package MNM, \emph{Journal of Statistical Software}, \bold{43}, 1-28.}
}
\author{Klaus Nordhausen}


\seealso{ \code{\link[ICSNP]{spatial.sign}}, \code{\link[SpatialNP]{spatial.rank}}, \code{\link[ICSNP]{HotellingsT2}}}

\examples{
X <- rmvt(150,diag(1,3))
g1 <- gl(3,50)
mv.Csample.test(X, g1)
mv.Csample.test(X, g1, score = "s")
mv.Csample.test(X, g1, score = "r")

Y <- rbind(rmvnorm(40,c(0,0,0)), rmvnorm(60,c(0,0,0.4)))
g2 <- factor(rep(1:2, c(40, 60)))
mv.Csample.test(Y, g2, score = "r")
mv.Csample.test(Y, g2, score = "r", method="p")
}
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}

