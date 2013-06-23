\docType{methods}
\name{open_well_response}
\alias{open_well_response}
\alias{open_well_response.default}
\title{Calculate the pressure/strain response spectrum for given formation properties (OPEN well)}
\usage{
  open_well_response(omega, T., S., Vw., Rs., Ku., B.,
    Avs. = 1, Aw. = 1, rho. = 1000, Kf. = 2.2e+09,
    grav. = 9.81, freq.units = NULL,
    model = c("liu", "cooper"))

  \method{open_well_response}{default} (omega, T., S., Vw.,
    Rs., Ku., B., Avs. = 1, Aw. = 1, rho. = 1000,
    Kf. = 2.2e+09, grav. = 9.81, freq.units = NULL,
    model = c("liu", "cooper"))
}
\arguments{
  \item{omega}{frequency, (see freq.units)}

  \item{T.}{effective aquifer transmissivity \eqn{[m^2/s]}}

  \item{S.}{well storativity, \eqn{[unitless]}}

  \item{Vw.}{well volume, \eqn{[m^3]}}

  \item{Rs.}{radius of screened portion, \eqn{[m]}}

  \item{Ku.}{undrained bulk modulus, \eqn{[Pa]}}

  \item{B.}{Skempton's coefficient, \eqn{[unitless,
  bounded]}}

  \item{Avs.}{amplification factor for volumetric strain
  \eqn{E_{kk,obs}/E_{kk}}, \eqn{[]}}

  \item{Aw.}{amplification factor of well volume change for
  \eqn{E_{kk}}, \eqn{[]}}

  \item{rho.}{fluid density \eqn{[kg/m^3]}}

  \item{Kf.}{bulk modulus of fluid, \eqn{[Pa]}}

  \item{grav.}{local gravitational acceleration
  \eqn{[m/s^2]}}

  \item{freq.units}{set what the units of frequency (omega)
  are: \code{"rad_per_sec"} (default, \code{NULL}), or
  \code{"Hz"}}

  \item{model}{character; use the model of Liu et al
  (1989), or Cooper et al (1965)}
}
\value{
  Matrix with three columns: radial frequency, amplitude,
  and phase [\eqn{\omega}), \eqn{A_\alpha (\omega)},
  \eqn{\Phi_\alpha (\omega)}] where the units of
  \eqn{\omega} will be radians per second, \eqn{A_\alpha
  (\omega)} in meters per strain, and \eqn{\Phi_\alpha
  (\omega)} in radians.
}
\description{
  Assumed values are: \describe{ \item{\code{Avs.}}{1}
  \item{\code{Aw.}}{1} \item{\code{rho.}}{\eqn{1000
  [kg/m^3]}, the density of water}
  \item{\code{Kf.}}{\eqn{2.2e9 [Pa]}, the bulk modulus of
  water} \item{\code{grav.}}{\eqn{9.81 [m/s^2]}, average
  gravitational force on Earth} }
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{well_response}}
}

