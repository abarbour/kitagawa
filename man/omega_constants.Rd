\docType{methods}
\name{omega_constants}
\alias{omega_constants}
\alias{omega_constants.default}
\title{calculate constants that depend on \eqn{\omega}}
\usage{
  omega_constants(omega = 0, c.type = c("alpha"), ...)

  \method{omega_constants}{default} (omega = 0,
    c.type = c("alpha"), ...)
}
\arguments{
  \item{omega}{frequency, \eqn{[rad/sec]}}

  \item{c.type}{the constant to calculate}

  \item{...}{additional params passed to calculator, e.g.
  S., T., Rs.}
}
\value{
  values of the 'c.type' constant (currently only 'alpha',
  i.e. eq 12)

  \code{NULL}
}
\description{
  e.g. Kitagawa equation 12
}
\examples{
omega_constants() # 0, with warnings about S., T., Rs.
omega_constants(T.=1,S.=1,Rs.=1)  # 0, no warnings
omega_constants(1:10)  # sequence, with warnings about S., T., Rs.
omega_constants(1:10,T.=1,S.=1,Rs.=1) # sequence, no warnings
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

