\name{summary.mvPCA}
\alias{summary.mvPCA}
\alias{print.summary.mvPCA}
\title{Summary for an object of class mvPCA.}
\description{Summary method for an object of class mvPCA.}
\usage{
\method{summary}{mvPCA}(object, loadings = FALSE, cutoff = 0.1, ...)

\method{print}{summary.mvPCA}(x, digits = 3, loadings = x$print.loadings, 
               cutoff = x$cutoff, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object to type \code{mvPCA}.}
  \item{loadings}{logical. Should the loadings be returned. }
  \item{cutoff}{ numeric. Loadings below this cutoff in absolute value are
          shown as blank in the output. }
  \item{x}{ an object of class "summary.mvPCA". }
  \item{digits}{the number of significant digits to be used in listing of
          loadings. }
  \item{\dots}{arguments to be passed to or from other methods. }
}

\value{
  'object' with additional components 'cutoff' and 'print.loadings'.
}

\author{ Klaus Nordhausen}
\seealso{ \code{\link{mvPCA}}}
\examples{
data(iris)
IRIS <- iris[,1:4]
iris.pca <- mvPCA(IRIS, "sign", "i")
summary(iris.pca, loadings = TRUE)
}
\keyword{print}
\keyword{methods}
