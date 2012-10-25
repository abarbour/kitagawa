\name{kitplot}
\alias{kitplot}
\title{Quickly plot the amplitude and phase spectra}
\usage{
  kitplot(Resp., xlim. = c(-4, 0),
    ylims = list(amp = c(5, 7), phs = 180 * c(-1, 1)),
    prep.resp = TRUE)
}
\arguments{
  \item{Resp.}{the response information as from
  well_response (see \strong{Details})}

  \item{xlim.}{frequency limits (assume log10 Hz scale)}

  \item{ylims}{list of limits (assume log10 m/strain for
  amplitude), use \code{amp} and \code{phs} for variable
  assignments}

  \item{prep.resp}{boolean, should the units of
  \code{Resp.} be appropriately transformed ?}
}
\value{
  A \code{data.frame} with the (transformed) values that
  have been plotted; this will include \emph{only} the
  subset of data falling in the range of \code{xlim.}
  inclusively.
}
\description{
  Used to mimic plots of frequency response spectra, as in
  Kitagawa et al (2011), e.g. \strong{Figures 7--9}.
}
\details{
  This is primarily a diagnostic tool, and is thus not very
  flexible in its implementation.

  The input data are assumed to be structured as they would
  be out of \code{\link{well_response}}, specifically three
  vectors representing: \describe{ \item{Radial
  frequency}{\eqn{[radians/s]}}
  \item{Amplitude}{\eqn{[m/strain]}}
  \item{Phase}{\eqn{[radians]}} } These will be transformed
  to have units of \eqn{[\log10 Hz]}, \eqn{[\log10
  m/strain]}, and \eqn{[rad]} respectively, unless
  \code{prep.resp=FALSE}.
}
\examples{
# dummy example: get some lines on the figure
n <- 10
ones <- rep(1, n)
fakeResp <- data.frame(f=2*pi*10**seq(-4, 0, length.out=n), amp=1e6*ones, phs=.9*pi*ones)
kitplot(fakeResp)
# focus in on a certain range:
fakeResp.foc <- kitplot(fakeResp, xlim.=c(-3, -1), ylims=list(amp=c(5.5, 6.5), phs=180*c(0, 1)))
kitplot(fakeResp.foc, prep.resp=FALSE) # fakeResp.foc has already been transformed
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\references{
  Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi
  (2011), Frequency characteristics of the response of
  water pressure in a closed well to volumetric strain in
  the high-frequency domain, \emph{J. Geophys. Res.},
  \strong{116}, B08301, doi:10.1029/2010JB007794.

  \url{http://www.agu.org/pubs/crossref/2011/2010JB007794.shtml}
}
\seealso{
  \code{\link{well_response}}
}

