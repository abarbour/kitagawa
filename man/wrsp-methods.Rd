\docType{methods}
\name{wrsp-methods}
\alias{as.data.frame.wrsp}
\alias{data.frame.wrsp}
\alias{kitplot}
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
  \method{as.data.frame}{wrsp} (x, ...)

  \method{data.frame}{wrsp} (x, ...)

  \method{print}{wrsp} (x, n = 3, ...)

  \method{summary}{wrsp} (object, ...)

  \method{print}{summary.wrsp} (x, ...)

  \method{lines}{wrsp} (x, series = c("amp", "phs"), ...)

  \method{points}{wrsp} (x, series = c("amp", "phs"),
    pch = "+", ...)

  \method{plot}{wrsp} (x, xlims = c(-3, 1),
    ylims = list(amp = NULL, phs = 185 * c(-1, 1)),
    logamp = TRUE, ...)

  kitplot(x, ...)

  \method{kitplot}{wrsp} (x, xlims = c(-3, 1),
    ylims = list(amp = NULL, phs = 185 * c(-1, 1)),
    logamp = TRUE, ...)
}
\arguments{
  \item{x}{'wrsp' object}

  \item{object}{'wrsp' object}

  \item{n}{numeric; the number of \code{head} and
  \code{tail} to print}

  \item{series}{character; the series to plot (amplitude or
  phase)}

  \item{pch}{point character, as in \code{\link{par}}}

  \item{xlims}{limits for x-axis (applies to both amp and
  phs frames)}

  \item{ylims}{optional list of limits for y-axis (i.e.,
  \code{list(amp=c(..),phs=c(...))})}

  \item{logamp}{logical; should the amplitude be in log10
  units}

  \item{...}{optional arguments}
}
\description{
  An object with class 'wrsp' is a list containing the
  response information, and the mechanical, hydraulic, and
  material properties used to generate the response for a
  sealed well.
}
\details{
  The response information is a matrix with frequency,
  complex response [\eqn{\omega}, \eqn{Z_\alpha (\omega)}]
  where the units of \eqn{\omega} will be as they were
  input. The amplitude of \eqn{Z} is in meters per strain,
  and the phase is in radians.

  \code{\link{kitplot}} was previously a standalone
  function, but is now simply a reference to
  \code{\link{plot.wrsp}}.
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
plot(Mod(wellresp) ~ omega, Wdf) # amplitude
plot(Arg(wellresp) ~ omega, Wdf) # phase
#
# or use the builtin method plot.wrsp
plot(W)
# change limits:
plot(W, xlims=c(-1,1), ylims=list(amp=c(5,8), phs=185*c(-1,1)))
}
\author{
  A. J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{well_response}}

  \code{\link{kitagawa-package}}

  Other PlotUtilities: \code{\link{log10_ticks}},
  \code{\link{log2_ticks}}, \code{\link{log_ticks}},
  \code{\link{logticks}}
}

