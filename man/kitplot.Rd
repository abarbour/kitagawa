\name{kitplot}
\alias{kitplot}
\title{Plot response spectrum}
\usage{
  kitplot(Resp., xlim = c(-4, 0),
    ylims = list(amp = c(5, 7), phs = 180 * c(-1, 1)))
}
\arguments{
  \item{Resp.}{the response information
  matrix(f,Amp,Phi,ncol=3)}

  \item{xlim}{frequency limits (assume log10 scale)}

  \item{ylims}{list of limits (assume log10 scale for
  amplitude), use 'amp' and 'phs'}
}
\value{
  NULL
}
\description{
  used to plot frequency response spectra as in kitagawa,
  e.g. Figures 7--9
}
\examples{
n <- 10
ones <- rep(1,n)
fakeResp <- data.frame(f=2*pi*10**seq(-4,0,length.out=n), amp=1e6*ones, phs=.9*pi*ones)
kitplot(fakeResp)
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

