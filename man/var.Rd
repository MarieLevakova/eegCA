\name{var}
\alias{var}
\title{
Estimation for a VAR model (if r = 0)
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
var(X, dt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
Data matrix, variables in columns
}
  \item{dt}{
Time step
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{A matrix out
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{Row 1}{Estimate of mu}
  \item{Row 2 to p+1}{Estimate of Omega}
  \item{Row p+2}{Trace test statistics}
  \item{Row p+3}{Eigenvalues}
  \item{Rows p+4 to p+N+4}{Residuals}
%% ...
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
