\name{vecm}
\alias{vecm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Core C++ routine for estimation parameters of VEC model by Johansen procedure
}
\description{
  The core of the Johansen procedure, implemented in C++.
}
\usage{
  vecm(X, r, A, B, dt)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A matrix containing the time series to be analyzed in its columns.
  }
  \item{r}{Cointegration rank}
  \item{A}{A matrix restricting alpha, such that alpha = A phi, where phi is an arbitrary vector.
  }
  \item{B}{A matrix restricting beta, such that beta = B psi, where psi is an arbitrary vector.
  }
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{A matrix ´out´
  %%  ~Describe the value returned
  %%  If it is a LIST, use
  \item{rows 1 to r}{Estimate of alpha}
  \item{rows r+1 to 2*r}{Estimate of beta}
  \item{row 2*r+1}{Estimate of mu}
  \item{row 2*r+2}{Trace test statistics}
  \item{row 2*r+3}{Eigenvalues}
  \item{rows 2*r+4 to 2*r+4+N}{Residuals}
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
