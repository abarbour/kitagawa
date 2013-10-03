\name{well_response}
\alias{well_response}
\alias{well_response.default}
\title{Spectral response for a sealed well}
\usage{
  well_response(omega, T., S., Vw., Rs., Ku., B., Avs, Aw,
    rho, Kf, grav, freq.units = c("rad_per_sec", "Hz"),
    as.pressure = TRUE)

  \method{well_response}{default} (omega, T., S., Vw., Rs.,
    Ku., B., Avs, Aw, rho, Kf, grav,
    freq.units = c("rad_per_sec", "Hz"),
    as.pressure = TRUE)
}
\arguments{
  \item{omega}{frequency, (see \code{freq.units})}

  \item{T.}{effective aquifer transmissivity \eqn{[m^2/s]}}

  \item{S.}{well storativity, \eqn{[unitless]}}

  \item{Vw.}{well volume, \eqn{[m^3]}}

  \item{Rs.}{radius of screened portion, \eqn{[m]}}

  \item{Ku.}{undrained bulk modulus, \eqn{[Pa]}}

  \item{B.}{Skempton's coefficient, \eqn{[unitless,
  bounded]}}

  \item{Avs}{amplification factor for volumetric strain
  \eqn{E_{kk,obs}/E_{kk}}, \eqn{[]}}

  \item{Aw}{amplification factor of well volume change for
  \eqn{E_{kk}}, \eqn{[]}}

  \item{rho}{fluid density \eqn{[kg/m^3]}}

  \item{Kf}{bulk modulus of fluid, \eqn{[Pa]}}

  \item{grav}{local gravitational acceleration
  \eqn{[m/s^2]}}

  \item{freq.units}{set the units of \code{omega}}

  \item{as.pressure}{logical; should the response for water
  pressure? (default is water heights)}
}
\value{
  An object with class 'wrsp'
}
\description{
  This is the primary function to calculate the response
  for a sealed well.
}
\details{
  The response depends strongly on the physical properties
  given. Default values are assumed where reasonable--for
  instance, the pore-fluid is assumed to be water--but
  considerable care should be invested in the choice of
  parameters, unless the function is used in an
  optimization scheme.

  Assumed values are: \tabular{rlrl}{ \code{Avs} \tab 1
  \tab \tab amplification factor for volumetric strain\cr
  \code{Aw} \tab 1 \tab \tab amplification factor for water
  well\cr }

  The responses returned here are, effectively, the
  amplification of water levels in a well, relative to the
  aquifer strain; or \deqn{Z = \frac{z}{\epsilon} \equiv
  \frac{p}{\rho g \epsilon}} If \code{as.pressure=TRUE},
  then the responses are scaled by \code{rho*grav} so that
  they represent water pressure relative to aquifer strain:
  \deqn{Z = \frac{p}{\epsilon}}

  Not all parameters need to be given, but should be. For
  example, if either \code{rho} or \code{grav} are not
  specified, they are taken from \code{\link{constants}}.
  \emph{Parameters which do not end in \code{.} do not need
  to be specified (they may be excluded); if they are
  missing, warnings will be thrown.}
}
\examples{
#### dummy example
well_response(1:10, T.=1, S.=1, Vw.=1, Rs.=1, Ku.=1, B.=1)

#### a more physically realistic calculation:
# Physical params applicable for B084 borehole
# (see: http://pbo.unavco.org/station/overview/B084/ for details)
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

# Access plot.wrsp:
plot(Rsp)
}
\author{
  A. J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{wrsp-methods}} for a description of the class
  'wrsp' and its methods, \code{\link{sensing_volume}} to
  easily estimate the volume \code{Vw.}, and
  \code{\link{kitagawa-package}} for references and more
  background.

  Other WellResponseFunctions:
  \code{\link{open_well_response}},
  \code{\link{open_well_response.default}}
}

