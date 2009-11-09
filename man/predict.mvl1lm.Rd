\name{predict.mvl1lm}
\Rdversion{1.1}
\alias{predict.mvl1lm}

\title{
Predicted Values Based on a Model Fitted by mv.l1lm 
}
\description{
Predicted response values based on a model fitted by \code{mv.l1lm}. 
}
\usage{
\method{predict}{mvl1lm}(object, newdata, na.action = na.pass, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
an object of class \code{mvl1lm}.
}
  \item{newdata}{
 An optional data frame with the values of the explaining variables.  If omitted, the fitted values are used.
}
  \item{na.action}{
function determining what should be done with missing values
          in 'newdata'.
}
  \item{\dots}{
needed for other methods.
}
}

\author{Klaus Nordhausen, \email{klaus.nordhausen@uta.fi} }

\keyword{methods}

