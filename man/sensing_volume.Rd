\name{sensing_volume}
\alias{sensing_volume}
\title{Calculate volume of fluids in the sensing borehole}
\usage{
  sensing_volume(rad_case, len_case, rad_screen,
    len_screen)
}
\arguments{
  \item{rad_case}{radius of the casing \eqn{[m]}}

  \item{len_case}{length of the casing \eqn{[m]}}

  \item{rad_screen}{radius of the screened interval
  \eqn{[m]}}

  \item{len_screen}{length of the screened interval
  \eqn{[m]}}
}
\value{
  scalar
}
\description{
  calculate Kitagawa equation 2
}
\examples{
sensing_volume(1, 1, 1, 1)  # dummy example for now
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}

