\name{plotShape}
\alias{plotShape}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pairwise Scatterplot Matrix of Shape Matrices}
\description{
  Function for visual comparisons for up to three shape matrices.
}
\usage{
plotShape(est1, est2 = NULL, est3 = NULL, X = NULL, alim = NULL, 
          color.ell = 2:4, color.points = grey(0.5), lty.ell = rep(1, 3), 
          pch.ell = rep(16, 3), pch.points = 1, level = 0.5, npoints = 100, 
          x.legend, y.legend, cex.legend = 1, pty = "s", gap = 1, 
          oma.bottom, labels, cex.labels = 2, main, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est1}{The shape matrix of interest. A list with the components \code{location}, \code{scatter} and \code{est.name}.}
  \item{est2}{An optional additional shape matrix. A list with the components \code{location}, \code{scatter} and \code{est.name}.}
  \item{est3}{An optional additional shape matrix. A list with the components \code{location}, \code{scatter} and \code{est.name}.}
  \item{X}{a numeric data frame or matrix. Optional data points on which the estimates could have been based.}
  \item{alim}{can be \code{NULL}, \code{both} or \code{ellipses}. Specifies when the plotting regions are computed if only the size of the ellipses are to be considered 
  or also the range of the data points.
            If \code{NULL} it chooses \code{both} if \code{X} is provided and otherwise \code{ellipses}.}
  \item{color.ell}{vector of length 3 that gives the colors for the corresponding estimates \code{'est1'}, \code{'est2'} and \code{'est3'}.}
  \item{color.points}{the color of the data points. }
  \item{lty.ell}{line types of the confidence ellipsoids.  }
  \item{pch.ell}{plotting symbols for the location estimates, the centers of the confidence ellipsoids. }
  \item{pch.points}{plotting symbol for the data points \code{X}.  }
  \item{level}{ The proportion of the data points that should be inside the ellipses. If there is no data the value for \code{t} in the function \code{ellipse}.}
  \item{npoints}{the number of points used to approximate each ellipsoid.}
  \item{x.legend}{vertical position of the legend. By default tries to find for 2 to 4 dimensional data a good location. If \code{NULL} no legend is drawn.}
  \item{y.legend}{horicontal position of the legend. By default tries to find for 2 to 4 dimensional data a good location. If \code{NULL} no legend is drawn.}
  \item{cex.legend}{\code{cex} for the legend.}
  \item{pty}{\code{pty} value for the individual plots of the scatter matrix. Default is \code{"s"}.}
  \item{gap}{distance between subplots, in margin lines.}
  \item{oma.bottom}{\code{oma} value of the bottom.}
  \item{labels}{optional labels for the diagonals. }
  \item{cex.labels}{\code{cex} for the labels. Default is 2. }
  \item{main}{ optional title of the plot. }
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  All scatter matrices are standardized to have determinant 1. If \code{X} is given, the Mahalanobis distances based on the location and shape estimates are computed, and \code{t}
  in the function \code{\link[ellipse]{ellipse} } is the level quantile of the Mahalanobis distances. If no \code{X} is provided \code{t} equals \code{level}.
  
  The location of the legend is currently problematic and it is
  recommended that the user should provide the coordinates for the legend.
}
\value{
  A scatter plot matrix.
}
\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }

\seealso{ \code{\link[ellipse]{ellipse}},  \code{\link{plotMvloc}}}

\examples{
X <- rmvt(100, diag(3), df=3)

EST1 <- list(location=colMeans(X), scatter=cov(X), est.name="COV")
HR <- HR.Mest(X)
EST2 <- list(location=HR$center, scatter=HR$scatter, est.name="Tyler")
plotShape(EST1,EST2, X=X)
}
\keyword{hplot }
