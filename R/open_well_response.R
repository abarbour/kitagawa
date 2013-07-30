#' Spectral response for an open well
#'
#' This is the primary function for an open (exposed to air) well.
#' 
#' @details
#' As opposed to \code{\link{well_response}}, this
#' calculates the theoretical, complex
#' well response for an unsealed (open) well.
#' 
#' The response depends strongly on the physical properties
#' given. Default values are assumed where reasonable--for instance, 
#' the pore-fluid is assumed to be water--but considerable care 
#' should be invested in the choice of
#' parameters, unless the function is used in an optimization scheme.
#' 
#' @section Models:
#' 
#' Rojstaczer (1988):
#' 
#' Cooper et al (1965) and Liu et al (1989):
#' These models are expressed succinctly in Roeloffs (1996).
#'
#' @name open_well_response
#' @export
#' 
#' @param omega  numeric; frequency,  (see \code{freq.units})
#' @param T. numeric; effective aquifer transmissivity \eqn{[m^2/s]}
#' @param S. numeric; well storativity,  \eqn{[unitless]}
#' @param z. numeric; the
#' @param freq.units charachter; setup the units of \code{omega}
#' @param model  character; setup the response model from 
#'    Rojstaczer (1988),
#'    Liu et al (1989), or 
#'    Cooper et al (1965)
#'
# @return Matrix with three columns: frequency, amplitude, and phase 
# [\eqn{\omega}), \eqn{A_\alpha (\omega)}, \eqn{\Phi_\alpha (\omega)}]
# where the units of \eqn{\omega} are as they were input,
# \eqn{A_\alpha (\omega)} in meters per strain, 
# and \eqn{\Phi_\alpha (\omega)} in radians.
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @seealso \code{\link{well_response}}, and
#' \code{\link{kitagawa-package}} for references and more background.
#'
open_well_response <- function(omega, T., S., z.,
                               freq.units=c("rad_per_sec","Hz"),
                               model=c("roj","liu","cooper")) UseMethod("open_well_response")
#' @rdname open_well_response
#' @method open_well_response default
#' @S3method open_well_response default
open_well_response.default <- function(omega, T., S., z.,
                                       freq.units=c("rad_per_sec","Hz"),
                                       model=c("rojstaczer","liu","cooper")){
  model <- match.arg(model)
  # Enforce units of omega to be radians/sec
  fc <- switch(match.arg(freq.units), rad_per_sec=1, Hz=2*pi)
  omega <- fc*omega
  #
  if (model=="rojstaczer"){
    Dtau. <- omega_constants(omega, c.type="diffusivity_time", S.=S., T.=T.)
    # Dtau == sqrt( omega / 2 / D.)
    # sqrt Qp <- z. * sqrt(omega / 2 / D) = z * Dtau
    sQp <- z. * Dtau.
    exptau <- exp(-sQp)
    A. <- exptau*cos(sQp) - 1
    B. <- -exptau*sin(sQp)
    wellresp <- complex(real=A., imaginary=B.)
  } else if (model=="liu" | model=="cooper"){
    
    # constants alpha and Kel
    Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
    Kmat <- alpha_constants(Alpha., c.type="Kel")
    Kel0. <- Kmat[,2] # zero order Kelvin function (complex)
    kei. <- Im(Kel0.)
    ker. <- Re(Kel0.)
    omegsq <- omega**2 # used a few times
    #
    # TODO
    # d is the aquifer thickness
    # Hw is the height of the water column above the upper limit of the aquifer
    Hw. <- d. <- 1
    #
    if (model=="liu"){
      .NotYetImplemented()
      #
      # from Liu et al (1989), expressed in Roeloffs 1996 eq 22
      #
      U. <- (d./T.)*Kel0.
      gamma <- sqrt(complex(imaginary=2*omega/(Rs.**2 * grav. * U.)))
      expgam <- exp(-1*gamma*d.)
      exp2gam <- exp(-2*gamma*d.)
      A. <- -1*omegsq/grav. * (Hw. + (1-expgam)/(1+expgam)/gamma)
      B. <- complex(imaginary=omega*U.*Rs.**2 * gamma*expgam/(1-exp2gam))
      wellresp <- 1/(A. - B. + 1)
      #
    } else if (model=="cooper"){
      .NotYetImplemented()
      #
      # from Cooper et al (1965) eq 28, expressed in Roeloffs 1996 eq 20
      #
      # the effective height of the water column in the well
      He. <- Hw. + 3*d./8 
      cT. <- omega * Rs.**2 / 2 / T.
      A. <- 1 - cT. * kei. + omegsq * He. / grav.
      B. <- cT. * ker.
      # the amplification of water level in the well relative to 
      # pressure head in the aquifer
      wellresp <- sqrt(A.**2 + B.**2)  # |w/h| or Mod(w/h)
      # eq 10
      # h/e_kk
      # h/e_areal ?
      # wellresp <- wellresp * -2 * mu * B. * (1+nu_u) /(rho. * grav * 3 * (1-2 * nu_u)) # | w/h * h/e_kk|
      omega_peak <- sqrt(grav./He.) # eq 21
    }
  }
  omega <- omega/fc
  toret <- cbind(omega, wellresp)
  return(toret)
}