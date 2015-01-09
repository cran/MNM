\name{mvPCA}
\alias{mvPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Principal Component Analysis }
\description{
  Principal component analysis based on different score functions
}
\usage{
mvPCA(X, score = "identity", estimate = "outer", 
      na.action = na.fail, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a numeric data frame or matrix with p columns. }
  \item{score}{score to be used. Can be either \code{"identity"}, \code{"sign"}, \code{"symmsign"} or \code{"rank"}.}
  \item{estimate}{can be \code{"outer"} or \code{"inner"}.}
  \item{na.action}{ a function which indicates what should happen when
    the data contain 'NA's. Default is to fail. }
  \item{\dots}{further arguments passed to or from other methods. }
}
\details{
  PCA as descriped in chapter 9 of the MNM book. Note that here ALL scatter matrices used are standardized to have trace(p).
  This function differs from most other PCA functions in R in that it does not center the data. The 'mvPCA' class has a print, summary, plot and predict
  method. 
}
\value{
  A list with class 'mvloc' containing the following components:
  \item{EigenV}{the standardized eigenvalues.}
  \item{loadings}{matrix with the corresponding loadings.}
  \item{scores}{matrix with the principal components.}
  \item{dname}{name of X.}
  \item{method}{Which shape matrix was used for the computation.}
  \item{n.obs}{number of observations used.}
  \item{p}{number of variables.}
}

\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } }
\author{Klaus Nordhausen}

\seealso{\code{\link{princomp}}, \code{\link{prcomp}}}

\examples{
data(iris)
IRIS <- iris[,1:4]
iris.pca <- mvPCA(IRIS, "sign", "i")
iris.pca
summary(iris.pca)
pairs(iris.pca$scores, col=iris[,5])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{nonparametric}
