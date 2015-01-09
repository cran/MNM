\name{predict.mvPCA}
\alias{predict.mvPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
 \title{ Prediction Method for a Principal Component Object of Type mvPCA} 
\description{Prediction method for class mvPCA.}
\usage{
\method{predict}{mvPCA}(object, newdata, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an object of class \code{mvloc}. }
  \item{newdata}{ New data with the same variables. If missing just the scores of object are returned.}
  \item{\dots}{ needed for other methods.}
}

\value{
 a matrix with the predicted principal components.
}
\author{Klaus Nordhausen}

\keyword{methods}
