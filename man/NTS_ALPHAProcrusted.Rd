\name{NTS_ALPHAProcrusted}
\alias{NTS_ALPHAProcrusted}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to estimate alpha
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
NTS_ALPHAProcrusted(Y, X, Z, ZBETA, r, Omega, P, beta, intercept = FALSE)
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
  \item{beta}{estimate of cointegrating vector}
  \item{intercept}{FALSE = do not include intercept, TRUE = include intercept in estimation short-run effects}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{ALPHA}{estimate of adjustment coefficients}
  \item{ALPHAstar}{estimate of transformed adjustment coefficients}
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
