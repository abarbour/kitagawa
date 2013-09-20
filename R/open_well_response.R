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
#' Rojstaczer (1988) is based on measurements of water level
#' and strain from volumetric or areal strainmeters.
#' 
#' Cooper et al (1965) and Liu et al (1989) are based
#' on measurements of water level and 
#' displacements from seismometers; these
#' models are expressed succinctly in Roeloffs (1996).
#'
#' @name open_well_response
#' @export
#' 
#' @param omega  numeric; frequency,  (see \code{freq.units})
#' @param T. numeric; effective aquifer transmissivity \eqn{[m^2/s]}
#' @param S. numeric; well storativity,  \eqn{[unitless]}
#' @param z. numeric; From Rojstaczer (1988): 'the depth from the water table'
#' @param Rs. numeric; the \emph{radius} of the open (screened) section
#' @param rho numeric; fluid density (assumed if missing)
#' @param grav numeric; the local gravitational acceleration (assumed if missing)
#' @param freq.units character; setup the units of \code{omega}
#' @param model  character; use the response model from 
#'    either
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
                               Rs.=(8/12)*(1200/3937),
                               rho, grav,
                               freq.units=c("rad_per_sec","Hz"),
                               model=c("rojstaczer","roj","liu","cooper")) UseMethod("open_well_response")
#' @rdname open_well_response
#' @method open_well_response default
#' @S3method open_well_response default
open_well_response.default <- function(omega, T., S., z.,
                                       Rs.=(8/12)*(1200/3937),
                                       rho, grav,
                                       freq.units=c("rad_per_sec","Hz"),
                                       model=c("rojstaczer","liu","cooper","hsieh")){
  model <- match.arg(model)
  # Enforce units of omega to be radians/sec
  fc <- switch(match.arg(freq.units), rad_per_sec=1, Hz=2*pi)
  omega <- fc*omega
  #
  const <- kitagawa::constants(FALSE)
  if (missing(rho)) rho <- const$water$density
  if (missing(grav)) grav <- const$gravity
  rhog <- rho*grav
  #
  Dtau. <- omega_constants(omega, c.type="diffusivity_time", S.=S., T.=T.)
  # Dtau == sqrt( omega / 2 / D.) == sqrt (omega * S / 2 / T)
  # sqrt Qp <- z. * sqrt(omega / 2 / D) = z * Dtau
  #
  if (model=="rojstaczer"){
    sQp <- z. * Dtau.
    exptau <- exp(-sQp)
    # scale by rhog for pressure over strain relative to static response
    A. <- rhog*(exptau*cos(sQp) - 1)
    B. <- -1*rhog*exptau*sin(sQp)
    wellresp <- complex(real=A., imaginary=B.)
  } else if (model %in% c("liu","cooper","hsieh")){
    # constants alpha and Kel
    Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
    KelPhiPsi <- alpha_constants(Alpha., c.type="Phi") #alpha, kel0, kel1, phi, psi
    Kel0. <- KelPhiPsi[,2] # zero order Kelvin function (complex)
    Phi. <- KelPhiPsi[,4] # zero order Kelvin function (complex)
    Psi. <- KelPhiPsi[,5] # zero order Kelvin function (complex)
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
    } else if (model=="hsieh"){
      #
      # R = rc**2 * omega / 2 T
      R. <- Rs.**2 * (Dtau.**2 / S.)
      U1. <- Psi. * ker.  +  Phi. * kei.
      U2. <- Phi. * ker.  +  Psi. * kei.
      A. <- 1  -  R. * U1.
      B. <- -1 * R. * U2.
      wellresp <- complex(real=A., imaginary=B.)
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
      #
      ##omega_peak <- sqrt(grav./He.) # eq 21
    }
  }
  omega <- omega/fc
  toret <- cbind(omega, wellresp)
  return(toret)
}