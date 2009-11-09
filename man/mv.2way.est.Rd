\name{mv.2way.est}
\alias{mv.2way.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Treatment Effect Estimates in the Randomized Complete Block Case}
\description{The treatment effect estimates for different
  score functions and their asymptotic covariance matrices in the
  randomized complete block case.}
\usage{
mv.2way.est(x, block, treatment, score = c("identity", "sign", "rank"),
             stand = c("outer", "inner"), 
             eps=1.0e-10, n.iter=1000, na.action = na.fail)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a numeric data frame or matrix.}
  \item{block}{a factor with at least two levels.}
  \item{treatment}{a factor with at least two levels.}
  \item{score}{the score to be used. Possible choices are
    \code{identity}, \code{sign} and \code{rank}.}
  \item{stand}{the standardization method used. Possible choices
    are \code{outer} and \code{inner}.}
  \item{eps}{convergence criterion.}
  \item{n.iter}{maximum number of iterations.}
  \item{na.action}{a function which indicates what should happen when
    the data contain 'NA's. Default is to fail.}
}
\details{This implements the treatment effect estimates described in
  chapter 12 of the MNM book.}

\value{
  A list of length (\code{treatment}-1) with class 'mvcloc'.
  The ith component of the previous list is a list with class 'mvloc'
  containing the following components:
  \item{location}{the treatment effect estimate when comparing the first
  treatment to the (i+1)th treatment.}
  \item{vcov}{the asymptotic covariance matrix of the treatment effect estimate.}
  \item{est.name}{name of the treatment effect estimate.}
  \item{dname}{name of data set for which the treatment effect estimate was computed.}
}

\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } }

\author{Jyrki \enc{Möttönen}{Mottonen} \email{jyrki.mottonen@helsinki.fi}}

\seealso{\code{\link{mv.2way.test}}, \code{\link{mv.1sample.est}}, \code{\link{mv.2sample.est}}}  

\examples{
data(beans)
est<-mv.2way.est(beans[,3:5],beans$Block,beans$Treatment,score="r",stand="i")
summary(est)
}
\keyword{multivariate}
\keyword{nonparametric}
