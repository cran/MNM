\name{spatial.sign2}
\Rdversion{1.1}
\alias{spatial.sign2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Spatial Signs
}
\description{
The function computes the spatial signs for a data set. This function differs from the function \code{\link[ICSNP]{spatial.sign}} in the way
how observations with small norms are treated. For details see below.
}
\usage{
spatial.sign2(X, center = TRUE, shape = TRUE, eps.S = 1e-05, 
              na.action = na.fail, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
a numeric data frame or matrix.
}
  \item{center}{
either a logical value or a numeric vector of length equal to
          the number of columns of 'X'. See below for more information. 
}
  \item{shape}{
either a logical value or a square numeric matrix with number of columns equal to
          the number of columns of 'X'. See below for more information.
}
  \item{eps.S}{
treshold value which defines which obersvations are considered to have a small norm.
}
  \item{na.action}{
a function which indicates what should happen when the data
          contain 'NA's.  Default is to fail.
}
  \item{\dots}{
arguments that can be passed on to functions used for the estimation of location
        and shape.
}
}
\details{
The spatial signs U of X with location \eqn{\mu}{mu} and shape V are given by transforming the data points
  \eqn{z_i = (x_{i}-\mu)V^{-\frac{1}{2}}}{z_i=(x_i-mu)V^{-0.5}} and then computing
  \deqn{u_{i}=\frac{z_i}{\| z_i \|}.}{u_i = z_i\\||z_i||.} 
  If a numeric value is given as 'center' and/or 'shape' these are used as \eqn{\mu}{mu} and/or V in the above formula. 
  If 'center' and/or 'shape' are 'TRUE' the values for \eqn{\mu}{mu} and/or V are estimated, if 'FALSE' the origin is used as the 
  value of \eqn{\mu}{mu} and/or the identity matrix as the value of V.
  
  When the norm \eqn{\| z_i \|}{||z_i||} is 0 then the spatial sign is set usually to 0 as for example in the function 
  \code{\link[ICSNP]{spatial.sign}}. Here however if  the spatial designs are defined as 
  \deqn{u_{i}=\frac{z_i}{\| z_i \|} I(\| z_i \| > eps.S) + \frac{z_i}{eps.S} I(\| z_i \| \leq eps.S).}{u_i=z_i\\||z_i|| I(||z_i||>eps.S)+ z_i\\eps.S  I(||z_i||<=eps.S).}
}
\value{
 a matrix with the spatial signs of the data as rows or the univariate signs as a px1 matrix. The
     centering vector and scaling matrix used are returned as attributes
     'center' and 'shape'.
}

\author{
Klaus Nordhausen}


\seealso{
\code{\link[ICSNP]{spatial.sign}}, \code{\link[ICSNP]{HR.Mest}}
}
\examples{
# comparing spatial.sign and spatial.sign2
data(pulmonary)
head(spatial.sign2(pulmonary, c(-0.1099999,-0.12,-4.3),FALSE))
head(spatial.sign(pulmonary, c(-0.1099999,-0.12,-4.3),FALSE))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ multivariate }
\keyword{nonparametric}
