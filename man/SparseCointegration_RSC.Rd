\name{SparseCointegration_RSC}
\alias{SparseCointegration_RSC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sparse cointegration function used in determine_rank function
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SparseCointegration_RSC(p, Y, X, Z, r, alpha.init = NULL, beta.init, max.iter = 25, conv = 10^-3, lambda.gamma = 0.001, lambda_beta = matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1), rho.glasso = 0.5, cutoff = 0.8, intercept = F, glmnetthresh = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{number of lagged differences}
  \item{Y}{Response Time Series}
  \item{X}{Time Series in Differences}
  \item{Z}{Time Series in Levels}
  \item{r}{cointegration rank}
  \item{alpha.init}{initial value for adjustment coefficients}
  \item{beta.init}{initial value for cointegrating vector}
  \item{max.iter}{maximum number of iterations}
  \item{conv}{convergence parameter}
  \item{lambda.gamma}{tuning paramter for short-run effects}
  \item{lambda_beta}{tuning parameter for cointegrating vector}
  \item{rho.glasso}{tuning parameter for inverse error covariance matrix}
  \item{cutoff}{cutoff value time series cross-validation approach}
  \item{glmnetthresh}{tolerance parameter glmnet function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{BETAhat}{estimate of cointegrating vectors}
  \item{ALPHAhat}{estimate of adjustment coefficients}
  \item{ZBETA}{estimate of short-run effects}
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
