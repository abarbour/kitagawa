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
#' If \code{as.pressure=TRUE}, then the responses are scaled by
#' \code{rho*grav} so they represent hydrostatic pressure; if
#' either \code{rho} or \code{grav} are not specified, they
#' are taken from \code{\link{constants}}.
#' 
#' Not all parameters need to be given, but should be.  
#' \emph{The parameters which do not end in \code{.} do
#' not need to be specified (they may be excluded); if missing
#' then warnings will be thrown.}
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
#' @param Rs. numeric; the \emph{radius} of the open (screened) section
#' @param rho numeric; fluid density (assumed if missing)
#' @param grav numeric; the local gravitational acceleration (assumed if missing)
#' @param z numeric; From Rojstaczer (1988): the depth from the water table (assumed if missing and if needed)
#' @param Hw numeric; height of water column above confined surface (assumed if missing and if needed)
#' @param Ta numeric; thickness of aquifer (assumed if missing and if needed)
#' @param freq.units character; setup the units of \code{omega}
#' @param model  character; use the response model from 
#' @param as.pressure logical; should the units of water height be returned as pressure?
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
#' @author A. J. Barbour <andy.barbour@@gmail.com>
#' 
#' @seealso \code{\link{well_response}}, and
#' \code{\link{kitagawa-package}} for references and more background.
#' @family WellResponseFunctions
#'
open_well_response <- function(omega, T., S.,
                               Rs.=(8/12)*(1200/3937),
                               rho, grav, z, Hw, Ta,
                               freq.units=c("rad_per_sec","Hz"),
                               model=c("rojstaczer","liu","cooper","hsieh"),
                               as.pressure=TRUE) UseMethod("open_well_response")
#' @rdname open_well_response
#' @method open_well_response default
#' @S3method open_well_response default
open_well_response.default <- function(omega, T., S., 
                                       Rs.=(8/12)*(1200/3937),
                                       rho, grav, z, Hw, Ta,
                                       freq.units=c("rad_per_sec","Hz"),
                                       model=c("rojstaczer","liu","cooper","hsieh"),
                                       as.pressure=TRUE){
  # Pick a model
  model <- match.arg(model)
  # Enforce units of omega to be radians/sec
  freq.units <- match.arg(freq.units)
  fc <- switch(freq.units, rad_per_sec=1, Hz=2*pi)
  omega <- fc*omega
  #
  # Setup constants
  const <- kitagawa::constants(FALSE)
  if (missing(rho)){
    rho <- const$water$density
  }
  if (missing(grav)){
    grav <- const$gravity
  }
  # scale by rhog for pressure over strain, relative to static response
  rhog <- ifelse(as.pressure, rho*grav, 1)
  #
  # Diffusiv time in unified framework
  Dtau. <- omega_constants(omega, c.type="diffusivity_time", S.=S., T.=T.)
  # e.g., Rojstaczer 1988 Eq 11:
  #   Qp is
  #     z^2 omega / 2 / D
  #   and we want sqrt(Q')
  #   Dtau is 
  #     sqrt( omega / 2 / D.) == sqrt (omega * S / 2 / T)
  #   so
  #     sqrt(Qp) == z * sqrt(omega / 2 / D) == z * Dtau
  #
  if (model=="rojstaczer"){
    # Rojstaczer 1988
    # Eq A3 - A4
    # z is the xxx
    if (missing(z)){
      z <- 1
      warning("Depth from the water table 'z' not given. using default")
    }
    sQp <- z * Dtau.
    exptau <- exp(-sQp)
    #
    A. <- rhog*(exptau*cos(sQp) - 1)
    B. <- -1*rhog*exptau*sin(sQp)
    #
    wellresp <- complex(real=A., imaginary=B.)
    #
  } else if (model %in% c("liu","cooper","hsieh")){
    # 
    # TODO: rhog scaling
    # 
    # Calc various constants needed
    #
    Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
    KelPhiPsi <- alpha_constants(Alpha., c.type="Phi") #alpha, kel0, kel1, phi, psi
    Kel0. <- KelPhiPsi[,2] # zero order Kelvin function (complex)
    Phi. <- KelPhiPsi[,4] 
    Psi. <- KelPhiPsi[,5]
    kei. <- Im(Kel0.)
    ker. <- Re(Kel0.)
    omegsq <- omega**2 # used a few times
    #
    # Hw is the height of the water column above the upper limit of the aquifer
    if (missing(Hw)){
      Hw <- 1
      warning("water column height 'Hw' not given. using default")
    }
    # Ta is the aquifer thickness
    if (missing(Ta)){
      Ta <- 1
      warning("aquifer thickness 'Ta' not given. using default")
    }
    .NotYetTested <- function() warning("this model has not yet been verified.")
    #
    if (model=="liu"){
      .NotYetTested()
      #
      # from Liu et al (1989), expressed in Roeloffs 1996 eq 22
      #
      U. <- (Ta/T.)*Kel0.
      #4: In complex(imaginary = 2 * omega/(Rs.^2 * grav * U.)) :
      #  imaginary parts discarded in coercion
      gamma <- sqrt(2*omega/(Rs.**2 * grav * U.))
      expgam <- exp(-1*gamma*Ta)
      exp2gam <- exp(-2*gamma*Ta)
      A. <- -1*omegsq/grav * (Hw + (1-expgam)/(1+expgam)/gamma)
      #B. <- complex(imaginary=)
      B. <- omega*U.*Rs.**2 * gamma*expgam/(1-exp2gam)
      #wellresp <- 1 / (A. - B. + 1)
      wellresp <- (A. - B. + 1)
      #
    } else if (model=="hsieh"){
      .NotYetTested()
      #
      # from Hsieh et al (1987), Eq X
      #
      # R = rc**2 * omega / 2 T
      R. <- Rs.**2 * (Dtau.**2 / S.)
      U1. <- Psi. * ker.  +  Phi. * kei.
      U2. <- Phi. * ker.  +  Psi. * kei.
      A. <- 1  -  R. * U1.
      B. <- -1 * R. * U2.
      wellresp <- complex(real=A., imaginary=B.)
      #
    } else if (model=="cooper"){
      .NotYetTested()
      #
      # from Cooper et al (1965) eq 28, expressed in Roeloffs 1996 eq 20
      #
      # the effective height of the water column in the well
      He. <- Hw + 3*Ta/8 
      cT. <- omega * Rs.**2 / 2 / T.
      A. <- 1 - cT. * kei. + omegsq * He. / grav
      B. <- cT. * ker.
      # the amplification of water level in the well relative to 
      # pressure head in the aquifer
      #wellresp <- sqrt(A.**2 + B.**2)  # |w/h| or Mod(w/h)
      wellresp <- complex(real=A., imaginary=B.)
      # eq 10
      # h/e_kk
      # h/e_areal ?
      # wellresp <- wellresp * -2 * mu * B. * (1+nu_u) /(rho. * grav * 3 * (1-2 * nu_u)) # | w/h * h/e_kk|
      #
      ##omega_peak <- sqrt(grav/He.) # eq 21
      #
    }
  }
  omega <- omega/fc
  toret <- list(
    Response=cbind(omega, wellresp)
  )
  return(toret)
}