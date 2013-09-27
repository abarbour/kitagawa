\docType{methods}
\name{owrsp-methods}
\alias{as.data.frame.owrsp}
\alias{data.frame.owrsp}
\alias{lines.owrsp}
\alias{owrsp}
\alias{owrsp-methods}
\alias{plot.owrsp}
\alias{points.owrsp}
\alias{print.owrsp}
\alias{print.summary.owrsp}
\alias{summary.owrsp}
\title{Generic methods for objects with class \code{'owrsp'}.}
\usage{
  \method{as.data.frame}{owrsp} (x, ...)

  \method{data.frame}{owrsp} (x, ...)

  \method{print}{owrsp} (x, n = 3, ...)

  \method{summary}{owrsp} (object, ...)

  \method{print}{summary.owrsp} (x, ...)

  \method{lines}{owrsp} (x, series = c("amp", "phs"), ...)

  \method{points}{owrsp} (x, series = c("amp", "phs"),
    pch = "+", ...)

  \method{plot}{owrsp} (x, xlims = c(-3, 1),
    ylims = list(amp = NULL, phs = 185 * c(-1, 1)), ...)
}
\arguments{
  \item{x}{'owrsp' object}

  \item{object}{'owrsp' object}

  \item{n}{numeric; the number of \code{head} and
  \code{tail} to print}

  \item{series}{character; the series to plot (amplitude or
  phase)}

  \item{pch}{point character, as in \code{\link{par}}}

  \item{xlims}{limits for x-axis (applies to both amp and
  phs frames)}

  \item{ylims}{optional list of limits for y-axis (i.e.,
  \code{list(amp=c(..),phs=c(...))})}

  \item{...}{optional arguments}
}
\description{
  An object with class 'owrsp' is a list containing the
  response information, and the mechanical, hydraulic, and
  material properties used to generate the response for an
  open well.
}
\details{
  The response information is a matrix with frequency,
  complex response [\eqn{\omega}, \eqn{Z_\alpha (\omega)}]
  where the units of \eqn{\omega} will be as they were
  input. The amplitude of \eqn{Z} is in meters per strain,
  and the phase is in radians.
}
\author{
  A. J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{open_well_response}}

  \code{\link{kitagawa-package}}
}

