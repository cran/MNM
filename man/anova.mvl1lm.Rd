\name{anova.mvl1lm}
\Rdversion{1.1}
\alias{anova.mvl1lm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Comparisons between Multivariate Linear Models
}
\description{
Comparisons between nested multivariate linear models fitted by \code{\link{mv.l1lm}}. The comparison can be based on score type of tests and Wald type of tests.
}
\usage{
\method{anova}{mvl1lm}(object, object2 = NULL, test = "Score", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{mvl1lm}. This gives the full model.
}
  \item{object2}{
an object of class \code{mvl1lm} or \code{NULL}. This gives the restricted (nested) model.
}
  \item{test}{
The test to be used. Options are \code{Score} and \code{Wald}. The score version is the default.
}
  \item{\dots}{
needed for other methods.
}
}
\details{
If only \code{object} is provided the function tests if all parameters equal zero. If \code{object} and \code{object2} are provided the function tests the
null hypothesis that the the restricted model (\code{object2}) is true.
For details see chapter 13 of the MNM book. Note that it is the users responsibility to make sure that the two models are nested and fitted on the same data.
For the regular L2 regression \code{\link{anova.mlm}} provides more options. 
}
\value{
A list with class 'anovamvl1lm' containing the following components:
  \item{models}{the model call(s) of \code{object} and \code{object2}.}
  \item{method}{type of the test used.}
  \item{statistic}{value of the test statistic.}
  \item{parameter}{degrees of freedom.}
  \item{p.value}{p-value of the test.}
}
  

\references{
\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer.}

\cite{Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: The Package MNM, \emph{Journal of Statistical Software}, \bold{43}, 1-28.} 
}
\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }

\examples{
# creating simple data

X <- cbind(rep(1,100),rmvnorm(100,c(0,0,0)) )
B <- matrix(c(4,1,1,0.5,-3,2,2,2),ncol=4, byrow=TRUE)
Y <- X \%*\% t(B)+ rmvnorm(100,c(0,0), diag(0.2,2))
DAT <- data.frame(x1=X[,2],x2=X[,3], x3=X[,4])

FullModel <- mv.l1lm(Y ~ x1 + x2 + x3, scores= "s", stand="i", data=DAT)
RestModel <- mv.l1lm(Y ~ x1, scores= "s", stand="i", data=DAT)

anova(FullModel)
anova(FullModel, RestModel)
anova(FullModel, RestModel, test="W")
}
\keyword{methods}
