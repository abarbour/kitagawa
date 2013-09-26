\name{open_well_response}
\alias{open_well_response}
\alias{open_well_response.default}
\title{Spectral response for an open well}
\usage{
  open_well_response(omega, T., S., z.,
    Rs. = (8/12) * (1200/3937), rho, grav,
    freq.units = c("rad_per_sec", "Hz"),
    model = c("rojstaczer", "roj", "liu", "cooper"),
    as.pressure = TRUE)

  \method{open_well_response}{default} (omega, T., S., z.,
    Rs. = (8/12) * (1200/3937), rho, grav,
    freq.units = c("rad_per_sec", "Hz"),
    model = c("rojstaczer", "liu", "cooper", "hsieh"),
    as.pressure = TRUE)
}
\arguments{
  \item{omega}{numeric; frequency, (see \code{freq.units})}

  \item{T.}{numeric; effective aquifer transmissivity
  \eqn{[m^2/s]}}

  \item{S.}{numeric; well storativity, \eqn{[unitless]}}

  \item{z.}{numeric; From Rojstaczer (1988): 'the depth
  from the water table'}

  \item{Rs.}{numeric; the \emph{radius} of the open
  (screened) section}

  \item{rho}{numeric; fluid density (assumed if missing)}

  \item{grav}{numeric; the local gravitational acceleration
  (assumed if missing)}

  \item{freq.units}{character; setup the units of
  \code{omega}}

  \item{model}{character; use the response model from}

  \item{as.pressure}{logical; should the units of water
  height be returned as pressure?  either Rojstaczer
  (1988), Liu et al (1989), or Cooper et al (1965)}
}
\description{
  This is the primary function for an open (exposed to air)
  well.
}
\details{
  As opposed to \code{\link{well_response}}, this
  calculates the theoretical, complex well response for an
  unsealed (open) well.

  The response depends strongly on the physical properties
  given. Default values are assumed where reasonable--for
  instance, the pore-fluid is assumed to be water--but
  considerable care should be invested in the choice of
  parameters, unless the function is used in an
  optimization scheme.

  If \code{as.pressure=TRUE}, then the responses are scaled
  by \code{rho*grav} so they represent hydrostatic
  pressure; if either \code{rho} or \code{grav} are not
  specified, they are taken from \code{\link{constants}}.
}
\section{Models}{
  Rojstaczer (1988) is based on measurements of water level
  and strain from volumetric or areal strainmeters.

  Cooper et al (1965) and Liu et al (1989) are based on
  measurements of water level and displacements from
  seismometers; these models are expressed succinctly in
  Roeloffs (1996).
}
\author{
  A. J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{well_response}}, and
  \code{\link{kitagawa-package}} for references and more
  background.

  Other WellResponseFunctions: \code{\link{well_response}},
  \code{\link{well_response.default}}
}

