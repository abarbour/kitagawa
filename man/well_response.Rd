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

  \item{freq.units}{set what the units of frequency (omega)
  are: "rad_per_sec" (default, NULL), or "Hz"}
}
\value{
  Matrix with three columns: \eqn{\omega}, \eqn{A_\alpha
  (\omega)}, \eqn{\Phi_\alpha (\omega)}

  where the units of \eqn{\omega} will be radians per
  second

  \code{NULL}
}
\description{
  calculate Kitagawa equation 17
}
\examples{
#### dummy example
well_response(1:10, T.=1, S.=1, Vw.=1, Rs.=1, Ku.=1, B.=1)

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
Volw <- sensing_volume(Rc, Lc, Rs, Ls) # m**3, ~= 1.8
#
Frqs <- 10**seq.int(from=-4,to=0,by=0.1) # log10-space
head(Rsp <- well_response(omega=Frqs, T.=1e-6, S.=1e-5,
Vw.=Volw, Rs.=Rs, Ku.=40e9, B.=0.2, freq.units="Hz"))
#
kitplot(Rsp)
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

