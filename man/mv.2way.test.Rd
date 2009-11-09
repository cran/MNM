\name{mv.2way.test}
\alias{mv.2way.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Randomized Complete Block Design.}
\description{Multivariate tests for testing the null hypothesis that there is no treatment effect in a randomized complete block design
  using different scores.}
\usage{
mv.2way.test(x, block, treatment, score = c("identity", "sign", 
             "rank"), stand = c("outer", "inner"), 
             method = c("approximation", "permutation"),
             n.simu = 1000, eps=1.0e-10, n.iter=10000, 
             na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric data frame or matrix of response variables.}
  \item{block}{a factor with at least two levels.}
  \item{treatment}{a factor with at least two levels.}
  \item{score}{the score to be used. Possible choices are
    \code{identity}, \code{sign} and \code{rank}.}
  \item{stand}{the standardization method used. Possible choices are
    \code{outer} and \code{inner}.}
  \item{method}{method for the computation of the p-value for the
    spatial sign and spatial rank tests. Possible choices are
    \code{approximation} and \code{permutation}.}
  \item{n.simu}{number of simulated permutations if \code{method="permutation"}.}
  \item{eps}{convergence criterion.}
  \item{n.iter}{maximum number of iterations.}
  \item{na.action}{a function which indicates what should happen when
    the data contain 'NA's. Default is to fail.}
}
\details{This implements the tests described in chapter 12 of the MNM book.}
\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the test statistic.}
  \item{parameter}{the degrees of freedom for the test statistic or the
    number of replications in the simulation.}
  \item{p.value}{the p-value for the test.}
  \item{null.value}{the specified null hypothesis value of the location.}
  \item{alternative}{a character string with the value 'two.sided'.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data set
    and of the grouping vector.}
}
\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } }

\author{Jyrki \enc{Möttönen}{Mottonen} \email{jyrki.mottonen@helsinki.fi}}

\seealso{\code{\link{mv.1sample.test}}, \code{\link{mv.Csample.test}} ,  \code{\link{mv.2way.est}}} 
\examples{

blocks <- gl(10, 5)
treatments <- factor(rep(1:5, 10))
X <- rmvnorm(n = 50, mean = c(1,2,3), sigma = diag(3))
mv.2way.test(X, blocks, treatments, score="r", stand="i", method="a")
}
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}
