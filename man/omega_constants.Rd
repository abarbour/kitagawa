\docType{methods}
\name{omega_constants}
\alias{.wc_alpha.default}
\alias{omega_constants}
\alias{omega_constants.default}
\title{calculate constants that depend on \eqn{\omega}}
\usage{
  omega_constants(omega = 0, c.type = c("alpha"), ...)

  \method{omega_constants}{default} (omega,
    c.type = c("alpha"), ...)

  \method{omega_constants}{default} (omega, S., T., Rs.)
}
\arguments{
  \item{omega}{frequency, \eqn{[rad/sec]}}

  \item{c.type}{the constant to calculate (currently only
  'alpha', i.e. eq 12)}

  \item{...}{additional params passed to calculator (should
  not be used)}
}
\value{
  Matrix with

  \code{NULL}

  \code{NULL}
}
\description{
  e.g. Kitagawa equation 12
}
\examples{
omega_constants(1:10)  # dummy example for now
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

