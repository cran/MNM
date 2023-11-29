\name{mv.l1lm}
\Rdversion{1.1}
\alias{mv.l1lm}

\title{
 Linear Regression Based on Identity, Spatial Sign or Spatial Rank Scores
}

\description{
This function fits a multivariate linear regression model based on identity, spatial sign or spatial rank scores. Both inner and outer standardization are possible.
}
\usage{
mv.l1lm(formula, scores = "identity", stand = "outer", 
        maxiter = 1000, eps = 1e-06, eps.S = 1e-06, 
        x = TRUE, y = TRUE, data, subset, na.action)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
an object of class "formula" (or one that can be coerced to
that class): a symbolic description of the model to be
fitted. The left part of the formula (the response) must be a  n x p matrix with at least two columns.
}
  \item{scores}{
score to be used. Can be either \code{"identity"}, \code{"sign"} or \code{"rank"}.
}
  \item{stand}{
can be \code{"outer"} or \code{"inner"}.
}
  \item{maxiter}{
maximum number of iterations. Used only for \code{score = "sign"} and  \code{score = "rank"}.
}
  \item{eps}{
convergence tolerance. Used only for \code{score = "sign"} or  \code{score = "rank"}.
}
  \item{eps.S}{
lower limit for the residual norms. Used only for \code{score = "sign"}  or  \code{score = "rank"} in the iteration procedure to avoid to divide by a zero norm.
}
  \item{x}{
logical. Indicating whether the design matrix 'x' returned from the model matrix should be stored. Default is \code{TRUE}. Might be needed for example 
in the \code{anova} function.  
}
  \item{y}{
logical. Indicating whether the response matrix 'y' should be stored. Default is \code{TRUE}. 
}
  \item{data}{
an optional data frame, list or environment (or object
          coercible by 'as.data.frame' to a data frame) containing the
          variables in the model.  If not found in 'data', the
          variables are taken from 'environment(formula)', typically
          the environment from which 'mv.l1lm' is called.
}
  \item{subset}{
an optional vector specifying a subset of observations to be
          used in the fitting process.
}
  \item{na.action}{
a function which indicates what should happen when the data
          contain 'NA's. 
}
}
\details{
The theory behind this function is described in detail in Chapter 13 of the MNM book. 

For regular multivariate L2-regression the function \code{\link{lm}} might be more efficient and offers more methods. Note however that
the results given by \code{\link{lm}} and \code{\link{mv.l1lm}} may differ slightly due to different divisors of the covariance matrix.

The algorithms for the sign and rank scores are still in an early phase and therefore any feedback is very welcome. 
For example if p+1 residuals are 0, then the algorithms may not return correct values. Note also that the computations for rank scores 
might be slow. 

Rank regression does not provide an estimate for the intercept parameter is not considered a parameter, a Hodges-Lehmann estimator of the residuals is then an estimate when an interecept term is in the formula.
For the one sample case however the function cannot be used for rank scores. We recommend that the regression function should not be used for the one or two sample case. There are
distinct functions designed for that purpose. Note furthermore that in the two sample case the covariance matrix returned from the regression function differs
slightly from the one returned by the function \code{\link{mv.2sample.est}} since there matrix A is computed in a different way. 

In general it is recommended to use the \code{data} argument and specify there the data frame that contains the variables and matrices. For having a matrix Y in a data frame 
for example the following methods work:
\describe{
\item{a)}{
\code{MyData <- data.frame(I(Y),...)}

or
}
\item{b)}{
\code{MyData <- data.frame(...)} \cr
\code{MyData$Y <- Y}
}}
Otherwise also the function \code{\link{cbind}} can be used on the left side of the formula to combine numeric vectors on the fly.
}
\value{

     \code{mv.l1ml} returns an object of 'class' \code{mvl1lm}.

     The functions \code{summary} is the best choice to view the results.  The
     generic accessor functions \code{coefficients}, \code{fitted}, \code{residuals}
     and \code{vcov} extract various useful features of
     the value returned by \code{mv.l1ml}.

     An object of class \code{mv.l1ml} is a list wich contains different information depending on the scores and standardization used.
     To see its content use the function \code{\link{str}}.
     
}


\references{
\cite{Oja, H. (2010), \emph{Multivariate Nonparametric Methods with R}, Springer.} 

\cite{Nordhausen, K. and Oja, H. (2011), Multivariate L1 Methods: The Package MNM, \emph{Journal of Statistical Software}, \bold{43}, 1-28.}
}
\author{Klaus Nordhausen}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lm}}, \code{\link{mv.1sample.est}}, \code{\link{mv.1sample.test}}, \code{\link{mv.2sample.est}}, \code{\link{mv.Csample.test}}
}
\examples{
# creating simple data

X <- cbind(rep(1,100),rmvnorm(100,c(0,0,0)) )
B <- matrix(c(4,1,1,0.5,-3,2,2,2),ncol=4, byrow=TRUE)
Y <- X \%*\% t(B)+ rmvnorm(100,c(0,0), diag(0.2,2))
DAT <- data.frame(x1=X[,2],x2=X[,3], x3=X[,4], Y=I(Y))

# true B
t(B)

# example using identity scores
test1 <- mv.l1lm(Y ~ x1 + x2 + x3, data=DAT)

print(test1)
summary(test1)
coef(test1)
vcov(test1)
head(fitted(test1))
head(residuals(test1))

# example using outer sign scores
test2 <- mv.l1lm(Y ~ x1 + x2 + x3, scores= "s", data=DAT)

print(test2)
summary(test2)
coef(test2)
vcov(test2)
head(fitted(test2))
head(residuals(test2))

# example using inner sign scores
test3 <- mv.l1lm(Y ~ x1 + x2 + x3, scores= "s", stand="i", 
data=DAT)

print(test3)
summary(test3)
coef(test3)
vcov(test3)
head(fitted(test3))
head(residuals(test3))

# example using outer rank scores
test4 <- mv.l1lm(Y ~ x1 + x2 + x3, scores= "r", stand="o", 
data=DAT)

print(test4)
summary(test4)
coef(test4)
vcov(test4)
head(fitted(test4))
head(residuals(test4))

# example using inner rank scores
test5 <- mv.l1lm(Y ~ x1 + x2 + x3, scores= "r", stand="i", 
data=DAT)

print(test5)
summary(test5)
coef(test5)
vcov(test5)
head(fitted(test5))
head(residuals(test5))

# prediction

newData <- data.frame(x1=c(1,-2),x2=c(0.5,0.7), x3=c(-1,-1))
newData
predict(test1,newData)
predict(test2,newData)
predict(test3,newData)
predict(test4,newData)
predict(test5,newData)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{multivariate}
\keyword{nonparametric}
\keyword{regression}
