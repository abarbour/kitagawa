\name{omega_norm}
\alias{omega_norm}
\title{Dimensionless frequency from diffusivity and depth}
\usage{
  omega_norm(omega, Diffusiv, z, invert = FALSE)
}
\arguments{
  \item{omega}{numeric; angular frequency}

  \item{Diffusiv}{numeric; hydraulic diffusivity}

  \item{z}{numeric; depth}

  \item{invert}{logical; should \code{omega} be taken as
  normalized frequency?}
}
\value{
  Dimensionless frequency unless \code{invert=TRUE}, which
  will return radial frequency.
}
\description{
  Dimensionless frequency from diffusivity and depth
}
\details{
  Dimensionless frequency \eqn{Q} is defined as
  \deqn{Q=\frac{z^2 \omega}{2 D}} where \eqn{z} is the well
  depth, \eqn{\omega} is the angular frequency and \eqn{D}
  is the hydraulic diffusivity
}
\seealso{
  \code{\link{open_well_response}},
  \code{\link{kitagawa-package}}

  Other utilities: \code{\link{.in0to1}},
  \code{\link{.nullchk}}, \code{\link{is.wrsp}},
  \code{\link{sensing_volume}}
}

