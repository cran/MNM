\name{beans}
\alias{beans}
\docType{data}
\title{ Randomized Block Experiment of Plots of Beans}
\description{
  Results of a randomized block experiment in the Cook Islands involving the effect of six different treatments on plots of beans infested by the serpentine leaf miner insect. 
}
\usage{data(beans)}

\format{
  A data frame with 24 observations on the following 5 variables.
  \describe{
    \item{\code{Block}}{a factor with levels \code{1} \code{2} \code{3} \code{4}.}
    \item{\code{Treatment}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}.}
    \item{\code{y1}}{a numeric vector. The number of miners per leaf.}
    \item{\code{y2}}{a numeric vector. The weight of beans per plot (in kg).}
    \item{\code{y3}}{a numeric vector. 1/sin(sqrt(p)), where p is the proportion of leaves infested with borer.}
  }
}

\details{
The value of variable \code{y3} in Block 4 for Treatment 2 is an estimate of a missing value.
}

\source{
  Data courtesy of Dr. R. Fullerton.
}
\references{
  \cite{Seber, G. A. F. (1998), \emph{Multivariate Observations}, London: Arnold.  }}
  
\examples{
data(beans)
plot(beans)
}
\keyword{datasets}
