\name{oscillator}
\alias{oscillator}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulation of a vector oscillator
}
\description{
  Simulates p oscillator process and return their xy-coordinates.
}
\usage{
  oscillator(N, dt, z0, alpha, beta, omega, freq, lvl, arma::mat S_phi, S_gam, model = c("std", "win"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{Number of simulations}
  \item{dt}{Discretization timestep}
  \item{z0}{Vector of initial values}
  \item{alpha}{Loading matrix}
  \item{beta}{Cointegration matrix}
  \item{omega}{Vector of phase attractors}
  \item{freq}{Vector of intrinsic frequency in the case of the degenerate gamma process, freq = mu.}
  \item{lvl}{Vector of levels for gamma process}
  \item{S_phi}{Covariance matrix of noise for phi}
  \item{S_gam}{Covariance matrix of noise for gamma}
  \item{model}{String specifying the model; "std" for the degenerate amplitude process, "win" for a Winfree oscillator.}
}
\details{
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  Matrix, first column contains times, the remaining ones are cartesian coordinates of the oscillators.
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
