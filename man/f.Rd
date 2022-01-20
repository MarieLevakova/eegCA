\name{f}
\alias{f}
\title{Cointegration relationship vector for phi-process}
\usage{
f(phi, gam, alpha, beta, omega = 0)
}
\description{
Function returns cointegration relationship vector for phi-process in a linearly coupled oscillator in the form

f(phi, gamma) = alpha*beta' (phi - omega) + gamma
}
\arguments{
  \item{phi}{Vector of values of phi process}
  \item{gam}{Vector of values of gamma process}
  \item{alpha}{Loading matrix alpha}
  \item{beta}{Cointegration matrix beta}
  \item{omega}{Vector of attractor states}
}
\examples{
hello()
}
