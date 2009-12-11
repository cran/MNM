\name{plotMvloc}
\alias{plotMvloc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Function to Plot Multivariate Location Estimates and Their Confidence Ellipsoids. }
\description{
  Visually presents and compares different multivariate location estimates and their confidence ellipsoids.
}
\usage{
plotMvloc(est1, est2 = NULL, est3 = NULL, X = NULL, alim = NULL, 
          color.ell = 2:4, color.points = grey(0.5), 
          lty.ell = rep(1, 3), pch.ell = rep(16, 3), 
          lwd.ell = rep(1, 3), cex.ell = rep(1, 3), 
          pch.points = 1, level = 0.95, npoints = 100, 
          x.legend, y.legend, cex.legend = 1, pty = "s", gap = 1, 
          oma.bottom, labels, cex.labels = 2, main, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{est1}{The location of interest. A list with the components \code{location}, \code{vcov} and \code{est.name}, for example an object of class 'mvloc'.}
  \item{est2}{An optional additional location estimate. A list with the components \code{location}, \code{vcov} and \code{est.name}, for example an object of class 'mvloc'.}
  \item{est3}{An optional additional location estimate. A list with the components \code{location}, \code{vcov} and \code{est.name}, for example an object of class 'mvloc'.}
  \item{X}{a numeric data frame or matrix. Optional data points on which the estimates could have been based.}
  \item{alim}{can be \code{NULL}, \code{both} or \code{ellipses}. Specifies whether the plotting regions are are based on the confidence ellipsoids only or also the range of the data points.
            If \code{NULL} it chooses \code{both} if \code{X} is provided and otherwise \code{ellipses}.}
  \item{color.ell}{vector of length 3 that gives the colors for the corresponding estimates \code{'est1'}, \code{'est2'} and \code{'est3'}.}
  \item{color.points}{the color of the data points. }
  \item{lty.ell}{line types of the confidence ellipsoids.  }
  \item{pch.ell}{plotting symbols for the location estimates, the centers of the confidence ellipsoids. }
  \item{lwd.ell}{line width values of the confidence ellipsoids. }
  \item{cex.ell}{cex values for the location estimates, the centers of the confidence ellipsoids. }
  \item{pch.points}{plotting symbol for the data points \code{X}.  }
  \item{level}{ the level of the confidence ellipsoids.}
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
  The figure can be used to compare different multivariate location estimates. The location of the legend is currently problematic and it is
  recommended that the user should provide the coordinates for the legend.
}
\value{
   A scatterplot matrix.
}

\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }

\seealso{ \code{\link[ellipse]{ellipse}}, \code{\link{plotShape} }}

\examples{
X <- rmvt(50, diag(c(1, 2)), 3)

est1 <- mv.1sample.est(X)
est2 <- mv.1sample.est(X, "sign")
est3 <- mv.1sample.est(X, "rank", "inner")

plotMvloc(est1)
plotMvloc(est1, est2, est3, X, alim="b", lty.ell=1:3, pch.ell=14:16)
plotMvloc(est1, est2, est3, X, alim="e")
}
\keyword{hplot }

