\name{unwrapPhase}
\alias{unwrapPhase}
\title{Unwraps phase process of an oscillator.}
\usage{
  unwrapPhase(Z)
}
\description{
  Extracts phase processes from analytical signals and unwraps them.
}
\arguments{
  \item{Z}{Matrix with columns containing complex-valued analytical signals of p oscillators.
  The first columns should contain times describing the sampling of the process, 
  even columns should contain x-coordinates,
  odd columns (except for the first one) should contain y-coordinates.}
}
\examples{
  # inVec - input signal
  # inVec.time - the corresponding time stamps
  hVec <- hilbert(inVec)
  Z <- cbind(inVec.time, Re(hVec), Im(hVec))
  phi <- unwrapPhase(Z)
}
