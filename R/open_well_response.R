#' Spectral response for an open well
#'
#' This is the primary function for an open (exposed to air) well.
#' It calculates the theoretical, complex
#' well response from either Cooper et al (1965) or Liu et al (1989).
#' These models are expressed succinctly in Roeloffs (1996).
#' 
#' @details
#' The response depends strongly on the physical properties
#' given. Default values are assumed where reasonable--for instance, 
#' the pore-fluid is assumed to be water--but considerable care 
#' should be invested in the choice of
#' parameters, unless the function is used in an optimization scheme.
#' 
#' Assumed values are:
#' \tabular{rlrl}{
#' \code{Avs.}  \tab 1 \tab \tab amplification factor for volumetric strain\cr
#' \code{Aw.}   \tab 1 \tab \tab amplification factor for water well\cr
#' \code{rho.}  \tab \eqn{1000}  \tab \eqn{[kg/m^3]} \tab the density of water \cr
#' \code{Kf.}   \tab \eqn{2.2}   \tab \eqn{[GPa]}    \tab the bulk modulus of water\cr
#' \code{grav.} \tab \eqn{9.81}  \tab \eqn{[m/s^2]}  \tab average gravitational force on Earth\cr
#' }
#' 
#' \emph{Note that Skempton's coefficient, \code{B.}, is bounded inclusively
#' within \eqn{[0,1]}; an error is thrown if it's not.
#' }
#' 
#' @name open_well_response
#' @export
#' 
#' @param omega  frequency,  (see \code{freq.units})
#' @param T.     effective aquifer transmissivity \eqn{[m^2/s]}
#' @param S.     well storativity,  \eqn{[unitless]}
#' @param Vw.    well volume,   \eqn{[m^3]}
#' @param Rs.    radius of screened portion,  \eqn{[m]}
#' @param Ku.    undrained bulk modulus,  \eqn{[Pa]}
#' @param B.     Skempton's coefficient,  \eqn{[unitless, bounded]}
#' @param Avs.   amplification factor for volumetric strain \eqn{E_{kk,obs}/E_{kk}},  \eqn{[]}
#' @param Aw.    amplification factor of well volume change for \eqn{E_{kk}},  \eqn{[]}
#' @param rho.   fluid density \eqn{[kg/m^3]}
#' @param Kf.    bulk modulus of fluid,  \eqn{[Pa]}
#' @param grav.  local gravitational acceleration \eqn{[m/s^2]}
#' @param freq.units  set the units of \code{omega}
#' @param model  character; use the model of Liu et al (1989), or Cooper et al (1965)
#'
#' @return Matrix with three columns: radial frequency, amplitude, and phase 
#' [\eqn{\omega}), \eqn{A_\alpha (\omega)}, \eqn{\Phi_\alpha (\omega)}]
#' where the units of \eqn{\omega} will be radians per second,
#' \eqn{A_\alpha (\omega)} in meters per strain, 
#' and \eqn{\Phi_\alpha (\omega)} in radians.
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @references Cooper, H. H., Bredehoeft, J. D., Papadopulos, I. S., and Bennett, R. R. (1965),
#' The response of well-aquifer systems to seismic waves, 
#' \emph{J. Geophys. Res.}, \strong{70} (16):3915-3926
#' 
#' @references Liu, L.-B., Roeloffs, E., and Zheng, X.-Y. (1989),
#' Seismically Induced Water Level Fluctuations in the Wali Well, Beijing, China,
#' \emph{J. Geophys. Res.}, \strong{94} (B7):9453-9462
#' 
#' @references Roeloffs, E. (1996),
#' Poroelastic techniques in the study of earthquake-related hydrologic phenomena,
#' \emph{Advances in Geophysics}, \strong{37}:135-195, Elsevier.
#' 
#' @seealso \code{\link{well_response}}
#'
open_well_response <- function(omega, T., S., Vw., Rs., Ku., B., 
                               Avs.=1,
                               Aw.=1,
                               rho.=1000, 
                               Kf.=2.2e9,
                               grav.=9.81,
                               freq.units=c("rad_per_sec","Hz"),
                               model=c("liu","cooper")) UseMethod("open_well_response")
#' @rdname open_well_response
#' @docType methods
#' @method open_well_response default
#' @S3method open_well_response default
open_well_response.default <- function(omega, T., S., Vw., Rs., Ku., B.,
                                       Avs.=1,
                                       Aw.=1,
                                       rho.=1000, 
                                       Kf.=2.2e9,
                                       grav.=9.81,
                                       freq.units=c("rad_per_sec","Hz"),
                                       model=c("liu","cooper")){
  # constants alpha and Kel
  Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
  Kmat <- alpha_constants(Alpha., c.type="Kel")
  Kel0. <- Kmat[,2] # zero order Kelvin function (complex)
  kei. <- Im(Kel0.)
  ker. <- Re(Kel0.)
  omegsq <- omega**2 # used a few times
  #
  model <- match.arg(model)
  #
  # TODO
  # d is the aquifer thickness
  # Hw is the height of the water column above the upper limit of the aquifer
  Hw. <- d. <- 1
  #
  if (model=="liu"){
    #
    # from Liu et al (1989), expressed in Roeloffs 1996 eq 22
    #
    U. <- (d./T.)*Kel0.
    gamma <- sqrt(complex(imaginary=2*omega/(Rs.**2 * grav. * U.)))
    expgam <- exp(-1*gamma*d.)
    exp2gam <- exp(-2*gamma*d.)
    A. <- -1*omegsq/grav. * (Hw. + (1-expgam)/(1+expgam)/gamma)
    B. <- complex(imaginary=omega*U.*Rs.**2 * gamma*expgam/(1-exp2gam))
    wellres <- 1/(A. - B. + 1)
    #
  } else if (model=="cooper"){
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
  toret <- cbind(omega, wellresp)
  return(toret)
}