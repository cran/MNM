\name{mv.1sample.est}
\alias{mv.1sample.est}

\title{ Multivariate One Sample Location Estimates}
\description{
  Estimates the multivariate location for different score functions and their asymptotic covariance matrices in the one sample case.
}
\usage{
mv.1sample.est(X, score = "identity", stand = "outer", maxiter = 100, 
               eps = 1e-06, na.action = na.fail, ... )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a numeric data frame or matrix. }
  \item{score}{ the score to be used. Possible choices are \code{identity}, \code{sign} and \code{rank}. }
  \item{stand}{ the standardization method used. Possible choices are \code{outer} and \code{inner}.}
  \item{maxiter}{ maximum number of iterations. Used only for \code{score = "sign"} and  \code{score = "rank"}.}
  \item{eps}{ convergence tolerance. Used only for \code{score = "sign"} and  \code{score = "rank"}.}
  \item{\dots}{ arguments that can be passed on to functions used for the estimation of location. }
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail. }
}
\details{
  For identity scores the location estimate is the regular mean vector. For the spatial sign score it is the spatial median in the outer
  standardization case and the Hettmansperger-Randles estimate in the inner standardization case. The rank estimate is the spatial Hodges-Lehmann estimator,
  either regular (\code{stand = "outer"}) or affine equivariant (\code{stand = "inner"}).
  
  Computation with outer standardization is faster than with inner standardization and especially the rank version might be slow and memory consuming.
  
  For further details see chapters 3, 5, 6, 7 and 8 of the MNM book. 
}
\value{
  A list with class 'mvloc' containing the following components:
  \item{location}{the location estimate as a vector.}
  \item{vcov}{the asymptotic covariance matrix of the location estimate.}
  \item{est.name}{name of the location estimate.}
  \item{dname}{name of the data set.}
}

\references{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. 

\cite{Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: The Package MNM, \emph{Journal of Statistical Software}, \bold{43}, 1-28.}}
\author{Klaus Nordhausen}


\seealso{ \code{\link[ICSNP]{spatial.sign}}, \code{\link[SpatialNP]{spatial.signrank}}, \code{\link[ICSNP]{spatial.median}}, \code{\link[ICSNP]{HR.Mest}}}
\examples{
  set.seed(1)
  X <- rmvt(100, diag(c(1, 2, 0.5)), 3)
 
  est.Hot.X <- mv.1sample.est(X)
  est.SS.o.X <- mv.1sample.est(X,"s")
  est.SS.i.X <- mv.1sample.est(X,"s","i")
  est.SR.o.X <- mv.1sample.est(X,"r")
  est.SR.i.X <- mv.1sample.est(X,"r","i")
  
  est.SR.o.X
  summary(est.SR.o.X)
  
  # plotting
   
  plot(est.Hot.X, est.SS.i.X, est.SR.i.X, X)
  # or
  plot(est.Hot.X, est.SS.i.X, est.SR.i.X)
}
\keyword{multivariate}
\keyword{nonparametric}
