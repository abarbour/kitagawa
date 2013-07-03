\docType{methods}
\name{well_response}
\alias{well_response}
\alias{well_response.default}
\title{Spectral response for a sealed well}
\usage{
  well_response(omega, T., S., Vw., Rs., Ku., B., Avs. = 1,
    Aw. = 1, rho. = 1000, Kf. = 2.2e+09, grav. = 9.81,
    freq.units = c("rad_per_sec", "Hz"))

  \method{well_response}{default} (omega, T., S., Vw., Rs.,
    Ku., B., Avs. = 1, Aw. = 1, rho. = 1000, Kf. = 2.2e+09,
    grav. = 9.81, freq.units = c("rad_per_sec", "Hz"))
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

  \item{Avs.}{amplification factor for volumetric strain
  \eqn{E_{kk,obs}/E_{kk}}, \eqn{[]}}

  \item{Aw.}{amplification factor of well volume change for
  \eqn{E_{kk}}, \eqn{[]}}

  \item{rho.}{fluid density \eqn{[kg/m^3]}}

  \item{Kf.}{bulk modulus of fluid, \eqn{[Pa]}}

  \item{grav.}{local gravitational acceleration
  \eqn{[m/s^2]}}

  \item{freq.units}{set the units of \code{omega}}
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
  This is the primary function for a sealed well. It
  calculates the theoretical, complex well response, namely
  \strong{Equation 17} in Kitagawa et al (2011). The
  results, however, are expressed as amplitude and phase.
}
\details{
  The response depends strongly on the physical properties
  given. Default values are assumed where reasonable--for
  instance, the pore-fluid is assumed to be water--but
  considerable care should be invested in the choice of
  parameters, unless the function is used in an
  optimization scheme.

  Assumed values are: \tabular{rlrl}{ \code{Avs.} \tab 1
  \tab \tab amplification factor for volumetric strain\cr
  \code{Aw.} \tab 1 \tab \tab amplification factor for
  water well\cr \code{rho.} \tab \eqn{1000} \tab
  \eqn{[kg/m^3]} \tab the density of water \cr \code{Kf.}
  \tab \eqn{2.2} \tab \eqn{[GPa]} \tab the bulk modulus of
  water\cr \code{grav.} \tab \eqn{9.81} \tab \eqn{[m/s^2]}
  \tab average gravitational force on Earth\cr }

  \emph{Note that Skempton's coefficient, \code{B.}, is
  bounded inclusively within \eqn{[0,1]}; an error is
  thrown if it's not. }
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
#
kitplot(Rsp)
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\references{
  Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi
  (2011), Frequency characteristics of the response of
  water pressure in a closed well to volumetric strain in
  the high-frequency domain, \emph{J. Geophys. Res.},
  \strong{116}, B08301, doi:10.1029/2010JB007794
}
\seealso{
  \code{\link{sensing_volume}} to estimate the volume
  \code{Vw.}

  \code{\link{kitplot}} to plot the results.

  \code{\link{open_well_response}} for modeling an open
  well.
}

