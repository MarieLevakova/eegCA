\name{johansen}
\alias{johansen}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Johansen procedure
}
\description{
A wrapper for Cpp function.
}
\usage{
johansen(X, r = 1, A = NULL, H = NULL, dt = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A matrix containing the time series to be analyzed in its columns.}
  \item{r}{Cointegration rank. If r = 0 (equivalent to stationarity), the function estimates a model with only a constant. Otherwise a standard VEC model assuming cointegration rank r is estimated.}
  \item{A}{Restrictions on alpha. The columns of A span a subspace where columns of alpha can live.}
  \item{H}{Restrictions on beta. The columns of H span a subspace where columns of beta can live.}
  \item{dt}{Time resolution of the observations. To be used for continuous stochastic processes, otherwise should take the default value 1.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of elements
  \item{N}{Number of observations minus 1.}
  \item{p}{Number of processes.}
  \item{r}{Cointegration rank (assigned at the function call).}
  \item{alpha}{Estimate of alpha matrix}
  \item{beta}{Estimate of beta matrix}
  \item{Psi}{Estimate of deterministic terms and coefficients of lags. Currently the function implements only a constant.}
  \item{Omega}{Estimate of the covariance matrix of innovations.}
  \item{test}{Test statistics of the rank test (maximum-eigenvalue form) of null hypotheses H(r=0), ..., H(r="r"-1)}
  \item{lambda}{Estimated eigenvalues.}
  \item{A}{Matrix of restrictions on alpha.}
  \item{H}{Matrix of restrictions on beta.}
  \item{df}{Degrees of freedom.}
  \item{dt}{Time resolution of the (continuous) process.}
  \item{res}{Residuals.}
  \item{data}{Matrix of submitted observations of the process.}
}

\references{
%% ~put references to the literature/web site here ~
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
