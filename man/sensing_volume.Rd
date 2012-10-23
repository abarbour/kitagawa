\name{sensing_volume}
\alias{sensing_volume}
\title{Calculate volume of fluids in the sensing borehole}
\usage{
  sensing_volume(rad_case, len_case, rad_screen,
    len_screen)
}
\arguments{
  \item{rad_case}{radius of the casing \eqn{[m]}}

  \item{len_case}{length of the casing \eqn{[m]}}

  \item{rad_screen}{radius of the screened interval
  \eqn{[m]}}

  \item{len_screen}{length of the screened interval
  \eqn{[m]}}
}
\value{
  scalar, the sensing volume in units of \eqn{[m^3]}
}
\description{
  calculate Kitagawa equation 2
}
\examples{
#### dummy example
sensing_volume(1, 1, 1, 1)
#
#### a more physically realistic calculation:
# Physical params applicable for B084 borehole
# (see: http://pbo.unavco.org/station/overview/B084/)
#
Rc <- 0.0508   # m, radius of water-sensing (2in)
Lc <- 146.9    # m, length of grouted region (482ft)
Rs <- 3*Rc     # m, radius of screened region (6in)
Ls <- 9.14     # m, length of screened region (30ft)
#
# calculate the sensing volume for the given well parameters
sensing_volume(Rc, Lc, Rs, Ls) # m**3, ~= 1.8
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

