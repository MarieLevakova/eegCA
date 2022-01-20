\name{NTS_BETA}
\alias{NTS_BETA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to estimate beta and Omega
}
\description{
First step: Determine Beta conditional on Gamma, alpha and Omega using Lasso Regression
Second step: Determine Omega conditional on Gamma, alpha and beta using GLasso
}
\usage{
NTS_BETA(Y, X, Z, ZBETA, r, Omega, P, alpha, alphastar, lambda_beta = matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1), rho.glasso, cutoff, intercept = F, glmnetthresh = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Response Time Series}
  \item{X}{Time Series in Differences}
  \item{Z}{Time Series in Levels}
  \item{ZBETA}{estimate of short-run effects}
  \item{r}{cointegration rank}
  \item{Omega}{estimate of inverse error covariance matrix}
  \item{P}{transformation matrix P derived from Omega}
  \item{alpha}{estimate of adjustment coefficients}
  \item{alphastar}{estimate of transformed adjustment coefficients}
  \item{lambda_beta}{tuning parameter for cointegrating vector}
  \item{rho.glasso}{tuning parameter for the inverse error covariance matrix}
  \item{cutoff}{cutoff value time series cross-validation approach}
  \item{intercept}{FALSE = do not include intercept, TRUE = include intercept in estimation short-run effects}
  \item{glmnetthresh}{tolerance parameter glmnet function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{BETA}{estimate of cointegrating vectors}
  \item{OMEGA}{estimate of inverse covariance matrix}
}
\references{
Wilms, I., Croux, C. (2016). Forecasting using sparse cointegration. International Journal of Forecasting, 32(4), 1256-1267.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
