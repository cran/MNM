\name{mv.2sample.est}
\alias{mv.2sample.est}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multivariate Two Sample Shift Estimates }
\description{
 Estimates the multivariate shift for different score functions and their asymptotic covariance matrices in the two sample case.
}
\usage{
mv.2sample.est(X, g, score = "identity", stand = "outer", 
               maxiter = 100, eps = 1e-06, na.action = na.fail,
               ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{ a numeric data frame or matrix. }
  \item{g}{ a factor with two levels }
  \item{score}{ the score to be used. Possible choices are \code{identity}, \code{sign} and \code{rank}. }
  \item{stand}{ the standardization method used. Possible choices are \code{outer} and \code{inner}.}
  \item{maxiter}{ maximum number of iterations. Used only for \code{score = "sign"} and  \code{score = "rank"}.}
  \item{eps}{ convergence tolerance. Used only for \code{score = "sign"} and  \code{score = "rank"}.}
  \item{na.action}{a function which indicates what should happen when the data
          contain 'NA's. Default is to fail.}
  \item{\dots}{ arguments that can be passed on to functions used for the estimation of location. }
}
\details{
  This implements the location estimates and their asymptotic covariance matrices as described in chapter 11 of the MNM book.
  Note that the shift is the parameter for the defference between 'values of level 1 - values of level 2' where the levels are as defined in the factor \code{g}.
  
  For the general c sample location case the function \code{\link{mv.l1lm}} might be used.  
}
\value{
  A list with class 'mvloc' containing the following components:
  \item{location}{the location estimate as a vector.}
  \item{vcov}{the asymptotic covariance matrix of the location estimate.}
  \item{est.name}{name of the location estimate.}
  \item{dname}{name of data set for which the location was computed.}
}
\references{\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer. } 

\cite{Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: The Package MNM, \emph{Journal of Statistical Software}, \bold{43}, 1-28.}
} 

\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }


\seealso{ \code{\link[ICSNP]{spatial.sign}}, \code{\link[SpatialNP]{spatial.signrank}}}

\examples{
X1<- rmvnorm(50,c(0,0,0))
X2<- rmvnorm(70,c(1,1,2))
X<-rbind(X1,X2)
g<-factor(rep(1:2,c(50,70)))

  est.Hot.X <- mv.2sample.est(X, g)
  est.SS.o.X <- mv.2sample.est(X, g, "s")
  est.SS.i.X <- mv.2sample.est(X, g, "s", "i")
  est.SR.o.X <- mv.2sample.est(X, g, "r")
  est.SR.i.X <- mv.2sample.est(X, g, "r", "i")

  est.SS.o.X
  
  summary(est.SS.o.X)
  
  # plotting
   
  plotMvloc(est.Hot.X, est.SS.i.X, est.SR.i.X)

}
\keyword{multivariate}
\keyword{nonparametric}
