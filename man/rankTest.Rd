\name{rankTest}
\alias{rankTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Test of the cointegration rank
}
\description{
Evaluates p-values of the rank test based on bootstrapped values.
}
\usage{
rankTest(M)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{M}{Output of the function bootstrap
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
List of elements
  \item{r}{Cointegration rank determined by the rank test}
  \item{pVal}{P-values of testing hypothesis H(r = 0), ..., H(r = p-1).}
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
