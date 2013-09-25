\name{omega_constants}
\alias{omega_constants}
\alias{omega_constants.default}
\title{Calculate any constants that depend on angular frequency \eqn{\omega}}
\usage{
  omega_constants(omega = 0,
    c.type = c("alpha", "diffusivity_time"), ...)

  \method{omega_constants}{default} (omega = 0,
    c.type = c("alpha", "diffusivity_time"), ...)
}
\arguments{
  \item{omega}{frequency, \eqn{[rad/sec]}}

  \item{c.type}{the constant to calculate}

  \item{...}{additional params passed to calculator.  For
  \code{ctype=="alpha"}, use \code{S., T., Rs.}, and for
  \code{ctype=="diffusivity_time"}, use \code{XXX}.}
}
\value{
  Values of the constant repesented by \code{c.type} for
  the given parameters
}
\description{
  This function accesses the appropriate method to
  calculate the \eqn{\omega}-dependent constant associated
  with the choice of \code{c.type}.
}
\details{
  \emph{This function is not likely to be needed by the
  user.}

  The name is in reference to radial frequency
  \eqn{\omega}, which is defined as \deqn{\omega \equiv 2
  \pi / \tau} where \eqn{\tau} is the period of
  oscillation.
}
\section{\code{"alpha"}

 The parameter \eqn{\alpha} is defined as
 \deqn{\alpha \equiv r_w\sqrt{\omega S / T}}
 where \eqn{r_w} is the radius of the well,
 where \eqn{S} is the storativity, and \eqn{T} is
 transmissivity.  See the parameter \code{...} for details
 on how to pass these values.

 This definition is common to many papers on the topic.  For example,
 it corresponds to  \strong{Equation 12} in Kitagawa et al (2011).
Because the computation of \eqn{\alpha} depends also on physical
properties, additional parameters can be
passed through (e.g. the transmissivity).}{

}

\section{\code{"diffusivity_time"}}{

}

\section{Warnings}{
  \strong{In the case \code{c.type='alpha'}(the default),
  the parameters \code{S.}, \code{T.}, and \code{Rs.} must
  be passed; otherwise, values are assumed to ensure the
  calculation does not fail, and the results are
  essentially meaningless.} Warnings will be issued if any
  necessary parameters are missing, indicating default
  values \code{S.=T.=Rs.=1} were used; these are physically
  unrealistic.
}
\examples{
# alpha
omega_constants() # default is alpha, but will give warnings about S., T., Rs.
omega_constants(T.=1,S.=1,Rs.=1)  # 0, no warnings
omega_constants(1:10)  # sequence, with warnings about S., T., Rs.
omega_constants(1:10,T.=1,S.=1,Rs.=1) # sequence, no warnings
# diffusivity time
omega_constants(c.type="diffusivity_time", D.=1)  # 0, no warnings
omega_constants(c.type="diff", D.=1)  # 0, no warnings (arg matching)
omega_constants(c.type="diff")  # 0, warnings about S., T. because no D.
omega_constants(c.type="diff", S.=1)  # 0, warnings about T. because no D. or S.
}
\author{
  Andrew Barbour <andy.barbour@gmail.com>
}
\seealso{
  \code{\link{alpha_constants}},
  \code{\link{well_response}}, and
  \code{\link{kitagawa-package}} for references and more
  background.
}

