\name{pairs2}
\Rdversion{1.1}
\alias{pairs2}

\title{
Plotting two numeric matrices
}
\description{
The function plots each variable contained in the matrix 'x' against the all variables contained in matrix 'y'. The function is not very sophisticated and 
only used to consider the residuals in a multivariate regression.
}
\usage{
pairs2(x, y, mars = c(4, 4, 0.1, 0.1), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
a numeric matrix. Same number or rows as \code{y}.
}
  \item{y}{
a numeric matrix. Same number or rows as \code{x}.
}
  \item{mars}{
A numerical vector of the form 'c(bottom, left, top, right)'
          which gives the number of lines of margin to be specified on
          the four sides of the plot. The default is \code{c(4, 4, 0.1, 0.1)}.
}
  \item{\dots}{
Arguments to be passed to methods, such as graphical
          parameters (see \code{\link{par}}). Should not contain \code{xlab} and \code{ylab}.
}
}


\author{Klaus Nordhausen}



\examples{
X <- rmvnorm(50, c(0,0,1))
Y <- rmvnorm(50, c(20,1), matrix(rep(0.5,4),ncol=2))
colnames(X) <- LETTERS[1:3]
colnames(Y) <- letters[1:2]
pairs2(X,Y)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot }

