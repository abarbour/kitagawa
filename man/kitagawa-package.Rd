\docType{package}
\name{kitagawa-package}
\alias{kitagawa}
\alias{kitagawa-package}
\title{Calculate theoretical spectral response for high-frequency borehole strain and fluid-pressure measurements}
\description{
  This package provides tools to calculate the theoretical
  spectral response of a closed well undergoing harmonic,
  volumetric strains. The calculation is based on the model
  published by Kitagawa et al (2011).
}
\details{
  The main function to be used is
  \code{\link{well_response}}, which takes in arguments for
  well- and aquifer-parameters, and frequencies at which to
  calculate the response functions. It accesses the
  constants-calculation routines (i.e.
  \code{\link{alpha_constants},\link{omega_constants}})
  where appropriate; hence, the user need not worry about
  those functions.

  There are also two helper functions included: \describe{
  \item{\code{\link{sensing_volume}}}{ to compute the
  sensing volume of fluid, for the specified well
  dimensions.} \item{\code{\link{kitplot}}}{ to subset and
  plot the resulting response information.  It is meant to
  mimic the format of figures published in Kitagawa et al
  (2011).} }
}
\section{Scientific background}{
  The model is based upon the assumption that fluid flows
  radially through an homogeneous, isotropic, confined
  aquifer.

  The underlying principle is as follows.  When a harmonic
  wave induces strain in a confined aquifer (one having
  aquitards above and below it), fluid flows radially into,
  and out of a well penetrating the aquifer. The
  flow-induced drawdown, \eqn{s}, is governed by the
  following partial differential equation, expressed in
  radial coordinates(\eqn{r}): \deqn{ \frac{\partial^2
  s}{\partial r^2} + \frac{1}{r} \frac{\partial s}{
  \partial r} - \frac{S}{T}\frac{\partial s}{\partial t} =
  0 } where \eqn{S, T} are the aquifer storativity and
  transmissivity respectively.

  The solution to this PDE, with periodic discharge
  boundary conditions, gives the amplitude and phase
  response we wish to calculate; this was originally
  presented by Hsieh et al (1987), and subsequently adapted
  by Kitagawa et al (2011) for the case of a sealed well.
  The model is applicable to any quasi-static process
  involving harmonic, volumetric strain of an aquifer (e.g.
  passing Rayleigh waves, or changes in the Earth's tidal
  potential)

  In practice, however, the presence of permeable fractures
  can violate the assumption of isotropic permeability,
  which may substantially alter the response by introducing
  shear-strain coupling. But these complications are beyond
  the scope of this model.
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\references{
  Hsieh, P. A., J. D. Bredehoeft, and J. M. Farr (1987),
  Determination of aquifer transmissivity from Earth tide
  analysis, \emph{Water Resour. Res.}, \strong{23} (10),
  1824-1832, doi:10.1029/WR023i010p01824.

  \url{http://www.agu.org/pubs/crossref/1987/WR023i010p01824.shtml}

  Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi
  (2011), Frequency characteristics of the response of
  water pressure in a closed well to volumetric strain in
  the high-frequency domain, \emph{J. Geophys. Res.},
  \strong{116}, B08301, doi:10.1029/2010JB007794.

  \url{http://www.agu.org/pubs/crossref/2011/2010JB007794.shtml}
}
\seealso{
  \code{\link{well_response}},
  \code{\link{sensing_volume}}, \code{\link{kitplot}}
}

