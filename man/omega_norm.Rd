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
  \code{\link{omega_norm}} returns dimensionless frequency,
  unless \code{invert=TRUE} where it will assume
  \code{omega} is dimensionless frequency, and return
  radial frequency.
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
\author{
  A. J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{open_well_response}},
  \code{\link{kitagawa-package}}

  Other utilities: \code{\link{sensing_volume}}
}

