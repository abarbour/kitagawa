\name{kitplot}
\alias{kitplot}
\title{Quickly plot the amplitude and phase spectra of sealed-well response}
\usage{
  kitplot(object, ...)
}
\arguments{
  \item{object}{the object to plot, with class 'wrsp', for
  example}

  \item{...}{additional parameters send to the plot method;
  \code{\link{plot.wrsp}}, for example}
}
\description{
  Reproduce plots of frequency response spectra, as in
  Kitagawa et al (2011, Figs 7--9).
}
\details{
  \emph{This is primarily a diagnostic tool, and is thus
  not very flexible in its implementation.}
}
\author{
  A. J. Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{well_response}} and
  \code{\link{wrsp-methods}}

  Other PlotUtilities: \code{\link{log10_ticks}},
  \code{\link{log2_ticks}}, \code{\link{log_ticks}},
  \code{\link{logticks}}
}

