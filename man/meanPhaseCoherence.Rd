\name{meanPhaseCoherence}
\alias{meanPhaseCoherence}
\title{Mean phase coherence measure}
\usage{
  meanPhaseCoherence(phi1, phi2)
}
\description{
  Calculates the mean phase coherence measure of a pair of oscillators with phase processes phi1 and phi2

R(phi1, phi2) = sqrt(Re(mean(exp(1i*(phi1-phi2))))^2+Im(mean(exp(1i*(phi1-phi2))))^2)
}
\arguments{
  \item{phi1}{Phase process of oscillator 1}
  \item{phi2}{Phase process of oscillator 2}
}
\examples{
  hello()
}
