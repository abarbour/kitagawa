\docType{methods}
\name{well_response}
\alias{well_response}
\alias{well_response.default}
\title{Calculate the pressure/strain response spectrum for given formation properties}
\usage{
  well_response(omega, T., S., Vw., Rs., Ku., B., Avs. = 1,
    Aw. = 1, rho. = 1000, Kf. = 2.2e+09, grav. = 9.81,
    freq.units = NULL)

  \method{well_response}{default} (omega, T., S., Vw., Rs.,
    Ku., B., Avs. = 1, Aw. = 1, rho. = 1000, Kf. = 2.2e+09,
    grav. = 9.81, freq.units = NULL)
}
\arguments{
  \item{omega}{frequency, (see freq.units)}

  \item{T.}{effective aquifer transmissivity \eqn{[m^2/s]}}

  \item{S.}{well storativity, \eqn{[]}}

  \item{Vw.}{well volume, \eqn{[m^3]}}

  \item{Rs.}{radius of screened portion, \eqn{[m]}}

  \item{Ku.}{undrained bulk modulus, \eqn{[Pa]}}

  \item{B.}{Skempton's coefficient, \eqn{[]}}

  \item{Avs.}{amplification factor for volumetric strain
  \eqn{E_{kk,obs}/E_{kk}}, \eqn{[]}}

  \item{Aw.}{amplification factor of well volume change for
  \eqn{E_{kk}}, \eqn{[]}}

  \item{rho.}{fluid density \eqn{[kg/m^3]}}

  \item{Kf.}{bulk modulus of fluid, \eqn{[Pa]}}

  \item{grav.}{local gravitational acceleration
  \eqn{[m/s^2]}}

  \item{freq.units}{set units of omega: "rad_per_sec"
  (default, NULL), or "Hz"}
}
\value{
  Matrix with three columns: \eqn{\omega}, \eqn{A_\alpha
  (\omega)}, \eqn{\Phi_\alpha (\omega)}

  \code{NULL}
}
\description{
  calculate Kitagawa equation 17
}
\examples{
well_response(1:10, 1, 1, 1, 1, 1, 1)  # dummy example for now
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

