\name{g}
\alias{g}
\title{Deterministic part of dGamma}
\usage{
  g(gam, lvl, model = c("std", "win"))
}
\description{
  Function g returns the deterministic part of dGamma either for a degenerate model or Winfree oscillator.

Degenerate model:
  g(gamma) = 0

Winfree oscillator:
  g(gamma) = (lvl - gamma)*gamma^2
}
\arguments{
\item{gam}{Vector of values of gamma process}
\item{lvl}{Attractor of the gamma process}
\item{model}{String specifying the model; "std" for a degenerate model, "win" for Winfree oscillator}
}
\examples{
hello()
}
