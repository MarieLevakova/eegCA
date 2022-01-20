\name{LRtest}
\alias{LRtest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Likelihood ratio test of restrictions on alpha and beta
}
\description{
Evaluates p-values of the likelihood ratio test of restrictions on alpha and beta.
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
  \item{from}{Unrestricted model (output of the function johansen)}
  \item{to}{Restricted model (output of the function johansen)}
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
