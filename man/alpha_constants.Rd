\docType{methods}
\name{alpha_constants}
\alias{alpha_constants}
\alias{alpha_constants.default}
\title{calculate constants that depend on \eqn{\alpha}}
\usage{
  alpha_constants(alpha = 0,
    c.type = c("Phi", "Psi", "A", "Kel"), ...)

  \method{alpha_constants}{default} (alpha = 0,
    c.type = c("Phi", "Psi", "A", "Kel"), ...)
}
\arguments{
  \item{alpha}{the constant alpha}

  \item{c.type}{the constant to calculate}

  \item{...}{additional params passed to calculator (should
  not be used)}
}
\value{
  Matrix with

  \code{NULL}
}
\description{
  e.g. # Kitagawa equations 10,11,18,19
}
\examples{
alpha_constants(1:10)  # dummy example for now
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

