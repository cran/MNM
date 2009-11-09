\name{screeplot.mvPCA}
\alias{plot.mvPCA}
\alias{screeplot.mvPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plotting Method for a Principal Component Object of Type mvPCA }
\description{
  Creates a screeplot for an object of class mvPCA. Works analogously to a normal screeplot for a classical principal component analysis.
  Here however the y-axis gives the proportion of the variation explained by the components.
}
\usage{
\method{plot}{mvPCA}(x, main = deparse(substitute(x)), ...)
\method{screeplot}{mvPCA}(x, npcs = min(10, length(x$EigenV)), 
          type = c("barplot", "lines"), 
          main = deparse(substitute(x)), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ an object to type \code{mvPCA}. }
  \item{npcs}{ the number of components to be plotted. }
  \item{type}{ the type of plot. }
  \item{main}{ title of the plot. }
  \item{\dots}{ other graphical parameters passed to or from other methods. }
}

\value{
  A screeplot.
}

\author{ Klaus Nordhausen, \email{klaus.nordhausen@uta.fi}}

\seealso{ \code{\link{mvPCA}}}
\examples{
data(IRIS)
IRIS <- iris[,1:4]
iris.pca <- mvPCA(IRIS, "sign", "i")
plot(iris.pca, type="lines")
}
\keyword{hplot}
\keyword{methods}
