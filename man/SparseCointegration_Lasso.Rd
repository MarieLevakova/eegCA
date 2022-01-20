\name{SparseCointegration_Lasso}
\alias{SparseCointegration_Lasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Main function to perform sparse cointegration
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SparseCointegration_Lasso(p, Y, X, Z, r, alpha.init = NULL, beta.init = NULL, max.iter = 10, conv = 10^-2, lambda.gamma = matrix(seq(from=0.01,to=0.001,length=5),nrow=1), rho.glasso=seq(from=1,to=0.1,length=5), lambda_beta=matrix(seq(from=0.1,to=0.001,length=100),nrow=1), cutoff=0.8, glmnetthresh=1e-04)
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
\item{lambda.gamma}{tuning paramter short-run effects}
\item{lambda_beta}{tuning paramter cointegrating vector}
\item{rho.glasso}{tuning parameter inverse error covariance matrix}
\item{cutoff}{cutoff value time series cross-validation approach}
\item{glmnetthresh}{tolerance parameter glmnet function}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
\item{BETAhat}{estimate of cointegrating vectors}
\item{ALPHAhat}{estimate of adjustment coefficients}
\item{ZBETA}{estimate of short-run effects}
\item{OMEGA}{estimate of inverse covariance matrix}
}

\references{
%% ~put references to the literature/web site here ~
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
