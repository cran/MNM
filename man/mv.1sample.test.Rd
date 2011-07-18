\name{mv.1sample.test}
\alias{mv.1sample.test}

\title{Multivariate Location Tests }
\description{
  Tests for multivariate location using different score functions.
}
\usage{
mv.1sample.test(X, mu = 0, score = "identity", stand = "outer", 
                method = "approximation", n.simu = 1000, 
                na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric data frame or matrix. }
  \item{mu}{ the null hypothesis value. Default is the zero vector.}
  \item{score}{the score to be used. Possible choices are \code{identity}, \code{sign} and \code{rank}.}
  \item{stand}{ the standardization method used. Possible choices are \code{outer} and \code{inner}. }
  \item{method}{ method for the computation of the p-value for the spatial sign and spatial signed-rank tests.
                Possible choices are \code{approximation} and \code{signchange}. }
  \item{n.simu}{ number of simulated sign changes if  \code{method=signchange}.}
  \item{na.action}{ a function which indicates what should happen when the data
          contain 'NA's. Default is to fail. }
}
\details{
 The tests provided here are the Hotelling's T^2 test, the spatial sign test and the signed-rank test and their affine invariant versions in the one sample location case. 
 
 Note that for the identity score the provided test is not the traditional Hotelling's T^2 test because here the covariance matrix is computed wrt to the null value and not wrt to the sample mean.
 Use the function \code{\link[ICSNP]{HotellingsT2}} for the traditional version of Hotelling's T^2 test.
 Details about the tests can be found in the chapters 5-8 of the MNM book.
}
\value{
 A list with class 'htest' containing the following components:
  \item{statistic }{the value of the test statistic.}
  \item{parameter}{the degrees of freedom for the test statistic or the number of replications in the simulation.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified hypothesized value of the location.}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data set.}
}
\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. }

\cite{Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: The Package MNM, \emph{Journal of Statistical Software}, \bold{43}, 1-28.}
 } 

\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }

 
\seealso{ \code{\link[ICSNP]{HotellingsT2}},  \code{\link[SpatialNP]{sr.loc.test}} }
\examples{
library(mvtnorm)
X <- rmvt(100, diag(c(1, 2, 0.5)), 3)

mv.1sample.test(X,mu=c(0,0,0.5))
mv.1sample.test(X,score="s", stand="i")
mv.1sample.test(X,score="s", stand="i", method="s")
mv.1sample.test(X,score="r", stand="o")
mv.1sample.test(X,score="r", stand="i")
}
\keyword{htest}
\keyword{ multivariate}
\keyword{nonparametric}
