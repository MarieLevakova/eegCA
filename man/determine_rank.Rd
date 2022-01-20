\name{determine_rank}
\alias{determine_rank}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function to determine the cointegration rank using the rank selection criterion
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
determine_rank(Y, X, Z, r.init = NULL, p, max.iter.lasso = 3, conv.lasso = 10^-2, max.iter.r = 5, beta.init, alpha.init, rho.glasso = 0.1, lambda.gamma = matrix(seq(from = 0.1, to = 0.001, length = 10), nrow = 1), lambda_beta = matrix(seq(from = 2, to = 0.001, length = 100), nrow = 1), glmnetthresh = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{Response Time Series}
  \item{X}{Time Series in Differences}
  \item{Z}{Time Series in Levels}
  \item{r.init}{initial value of cointegration rank}
  \item{p}{number of lags to be included}
  \item{max.iter.lasso}{maximum number of iterations PML}
  \item{conv.lasso}{convergence parameter}
  \item{max.iter.r}{maximum number of iterations to compute cointegration rank}
  \item{beta.init}{initial value for beta}
  \item{alpha.init}{initial value for alpha}
  \item{rho.glasso}{tuning parameter for inverse covariance matrix}
  \item{lambda.gamma}{tuning parameter for GAMMA}
  \item{lambda.beta}{tuning parameter for BETA}
  \item{glmnetthresh}{tolerance parameter glmnet function}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{rhat}{estimated cointegration rank}
  \item{it.r}{number of iterations}
  \item{rhat_iterations}{estimate of cointegration rank in each iteration}
  \item{mu}{value of mu}
  \item{decomp}{eigenvalue decomposition}
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
