\docType{methods}
\name{kitplot.wrsp}
\alias{as.data.frame.wrsp}
\alias{data.frame.wrsp}
\alias{kitplot.wrsp}
\alias{lines.wrsp}
\alias{plot.wrsp}
\alias{points.wrsp}
\alias{print.summary.wrsp}
\alias{print.wrsp}
\alias{summary.wrsp}
\alias{wrsp}
\alias{wrsp-methods}
\title{Generic methods for objects with class \code{'wrsp'}.}
\usage{
  \method{kitplot}{wrsp} (object, ...)

  \method{as.data.frame}{wrsp} (object, ...)

  \method{data.frame}{wrsp} (object, ...)

  \method{print}{wrsp} (object, ...)

  \method{summary}{wrsp} (object, ...)

  \method{print}{summary.wrsp} (object, ...)

  \method{lines}{wrsp} (object, series = c("amp", "phs"),
    ...)

  \method{points}{wrsp} (object, series = c("amp", "phs"),
    pch = "+", ...)

  \method{plot}{wrsp} (object, xlims = c(-3, 1),
    ylims = list(amp = NULL, phs = 185 * c(-1, 1)), ...)
}
\arguments{
  \item{object}{wrsp object}

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
  An object with class 'wrsp' is a list containing the
  response information, and the mechanical, hydraulic, and
  material properties used to generate the response.
}
\details{
  The response information is a matrix with three columns:
  frequency, amplitude, and phase [\eqn{\omega},
  \eqn{A_\alpha (\omega)}, \eqn{\Phi_\alpha (\omega)}]
  where the units of \eqn{\omega} will be as they were
  input, \eqn{A_\alpha (\omega)} is in meters per strain,
  and \eqn{\Phi_\alpha (\omega)} is in radians.
}
\examples{
W <- well_response(1:10, T.=1, S.=1, Vw.=1, Rs.=1, Ku.=1, B.=1)
str(W)
print(W)
print(summary(W))
#
# Plot the response
plot(rnorm(10), xlim=c(-1,11), ylim=c(-2,2))
lines(W)
lines(W, "phs", col="red")
points(W)
points(W, "phs")
#
Wdf <- as.data.frame(W)
plot(Amp. ~ omega, Wdf)
plot(Phs. ~ omega, Wdf)
#
# or use the builtin method
plot(W)
}
\author{
  A. J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{well_response}}

  \code{\link{kitagawa-package}}
}

