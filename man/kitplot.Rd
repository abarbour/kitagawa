\name{kitplot}
\alias{kitplot}
\title{Plot response spectrum}
\usage{
  kitplot(Resp., xlim = c(-4, 0),
    ylims = list(amp = c(5, 7), phs = 180 * c(-1, 1)))
}
\arguments{
  \item{Resp.}{the response information as from
  well_response: matrix(f,Amp,Phi,ncol=3)}

  \item{xlim}{frequency limits (assume log10 Hz scale)}

  \item{ylims}{list of limits (assume log10 m/strain for
  amplitude), use 'amp' and 'phs'}
}
\value{
  NULL
}
\description{
  Used to mimic plots of frequency response spectra, as in
  Kitagawa et al 2011, e.g. Figures 7--9
}
\examples{
# dummy example: get some lines on the figure
n <- 10
ones <- rep(1,n)
fakeResp <- data.frame(f=2*pi*10**seq(-4,0,length.out=n), amp=1e6*ones, phs=.9*pi*ones)
kitplot(fakeResp)
# focus:
kitplot(fakeResp, xlim=c(-3,1), ylims=list(amp=c(5.5,6.5),phs=180*c(0,1)))
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

