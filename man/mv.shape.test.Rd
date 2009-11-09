\name{mv.shape.test}
\alias{mv.shape.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Test for Sphericity }
\description{
 Test for sphericity based on different score functions.
}
\usage{
mv.shape.test(X, score = "identity", location = "est",
              na.action = na.fail, ...)
}
\arguments{
  \item{X}{ a numeric data frame or matrix.  }
  \item{score}{the score to be used. Possible are \code{identity}, \code{sign}, and \code{symmsign}.}
  \item{location}{ specifies if the location should be estimated or taken to be the origin. Possible choices are \code{est} and \code{origin}.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail.}
  \item{\dots}{ arguments  passed on to other functions.}
}
\details{
 Note that here inner standardization is not logical. The rank score test is not implemented. Otherwise the tests are as described in chapter
 9 of the MNM book. 
 
 To test for other "shapes" than sperical, transform the data accordingly and then test for spericity.
}
\value{
  A list with class 'htest' containing the following components:
  \item{statistic }{the value of the test statistic.}
  \item{parameter}{the degrees of freedom for the test statistic or the number of replications in the simulation.}
  \item{p.value}{the p-value for the test.}
  \item{method}{a character string indicating what type of test was performed.}
  \item{data.name}{a character string giving the name of the data used.}
}
\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } }
\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }


\seealso{\code{\link[SpatialNP]{sr.sphere.test}}}

\examples{
X <- rmvt(150,diag(1,3))
mv.shape.test(X)
mv.shape.test(X,"sym")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{multivariate}
\keyword{nonparametric}

