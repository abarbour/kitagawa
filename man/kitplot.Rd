\name{kitplot}
\alias{kitplot}
\title{Plot response spectrum}
\usage{
  kitplot(Resp., xlim = c(-4, 0))
}
\arguments{
  \item{Resp.}{the response information
  matrix(f,Amp,Phi,ncol=3)}

  \item{xlim}{frequency limits (assume log10 scale)}
}
\value{
  NULL
}
\description{
  used to plot frequency response spectra as in kitagawa,
  e.g. Figures 7--9
}
\examples{
kitplot(data.frame(f=2*pi*10**seq(-4,0,length.out=10), amp=1e6*rep(1,10), phs=.9*pi*rep(1,10))) # bad example
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

