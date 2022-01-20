\name{NTS_GAMMAFIXEDLassoReg}
\alias{NTS_GAMMAFIXEDLassoReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to estimate GAMMA, using L1 PENALTY, fixed value of lmabda.gamma
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NTS_GAMMAFIXEDLassoReg(Y, X, Z, Pi, Omega, lambda.gamma = 0.001, p, cutoff = 0.8, intercept = F, glmnetthresh = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Response Time Series}
  \item{X}{Time Series in Differences}
  \item{Z}{Time Series in Levels}
  \item{Pi}{estimate of cointegration space}
  \item{Omega}{estimate of inverse error covariance matrix}
  \item{lambda.gamma}{tuning parameter for short-run effects}
  \item{p}{number of lagged differences}
  \item{cutoff}{cutoff value time series cross-validation approach}
  \item{intercept}{FALSE = do not include intercept, TRUE = include intercept}
  \item{glmnetthresh}{tolerance parameter glmnet function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{ZBETA}{estimate of short-run effects}
  \item{P}{transformation matrix P derived from Omega}
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
