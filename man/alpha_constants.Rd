\docType{methods}
\name{alpha_constants}
\alias{alpha_constants}
\alias{alpha_constants.default}
\title{calculate constants that depend on \eqn{\alpha}}
\usage{
  alpha_constants(alpha = 0,
    c.type = c("Phi", "Psi", "A", "Kel"))

  \method{alpha_constants}{default} (alpha = 0,
    c.type = c("Phi", "Psi", "A", "Kel"))
}
\arguments{
  \item{alpha}{the constant alpha}

  \item{c.type}{the constant to calculate}
}
\value{
  Matrix with

  \code{NULL}
}
\description{
  e.g. # Kitagawa equations 10,11,18,19
}
\examples{
alpha_constants(0)  # 0 - Phi
alpha_constants(1:10)  # sequence - Phi
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

