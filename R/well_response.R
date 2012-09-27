#' Calculate the pressure/strain response spectrum for given formation properties
#'
#' calculate Kitagawa equation 17
#'
#' @name well_response
#' @export
#' 
#' @param omega  frequency,  (see @param freq.units)
#' @param rho.   fluid density,  \eqn{[kg/m^3]]
#' @param grav.  local gravitational acceleration, 	\eqn{[m/s^2]}
#' @param T.     effective aquifer transmissivity,  \eqn{[m^2/s]}
#' @param S.     well storativity,  \eqn{[]}
#' @param Vw.    well volume,	 \eqn{[m^3]}
#' @param Aw.    amplification factor of well volume change for \eqn{E_{kk}},  \eqn{[]}
#' @param Avs.   amplification factor for volumetric strain \eqn{E_{kk,obs}/E_{kk}},  \eqn{[]}
#' @param Ku.    undrained bulk modulus,  \eqn{[Pa]}
#' @param Kf.    bulk modulus of fluid,  \eqn{[Pa]}
#' @param Rs.    radius of screened portion,	\eqn{[m]}
#' @param B.     Skempton's coefficient,  \eqn{[]}
#' @param freq.units  units of @param omega: "rad_per_sec" (default, NULL), or "Hz"
#'
#' @return Matrix with three columns: \eqn{\omega}, \eqn{A_\alpha (\omega)}, \eqn{\Phi_\alpha (\omega)}
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' well_response(1:10, 1, 1, 1, 1, 1, 1)  # dummy example for now
well_response <-
function(omega,
         T., S., Vw., Rs., Ku., B.,
         Avs.=1,
         Aw.=1,
         rho.=1000, 
         Kf.=2.2e9,
         grav.=9.81,
         freq.units=NULL) UseMethod("well_response")

#' @return \code{NULL}
#' @rdname well_response
#' @docType methods
#' @method well_response default
#' @S3method well_response default
well_response.default <-
  function(omega, T., S., Vw., Rs., Ku., B., Avs., Aw.,
           rho., Kf., grav., freq.units){
    #
    # calculate Kitagawa equation 17
    #
    fc <- switch(match.arg(freq.units, c("rad_per_sec","Hz")),
                 rad_per_sec=1,
                 Hz=2*pi
    )
    omega <- fc*omega
    # Alpha function
    Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
    # A1, and A2 (functions of Phi and Psi, calculated internally)
    Amat <- alpha_constants(Alpha., c.type="A")
    stopifnot(ncol(Amat)==7)
    rm(Alpha.)  # cleanup
    #  A1,2 are in Mod(A.[,6:7]) 
    A12 <- Mod(Amat[,6:7])  # is complex, but imag is zero, so == abs
    rm(Amat)    # cleanup
    A1 <- A12[,1]
    A2 <- A12[,2]
    rm(A12)     # cleanup
    #
    # prevent duplicate calculations
    rhog <- rho. * grav.
    #print(summary(omega))
    TVFRG <- 2 * pi * T. / omega / Vw. / rhog
    #
    # calculate amp and phase of response
    #
    tmpd. <- Ku. * B. / Aw. * TVFRG - A2
    rNum. <- tmpd. * tmpd. + A1 * A1
    rm(tmpd.)
    tmpd. <- Kf. * TVFRG  -  A2
    rDen. <- tmpd. * tmpd. + A1 * A1
    rm(tmpd.)
    # amplitude, Kitagawa equation 20
    Amp. <- Kf. * Aw. / Avs. / rhog * sqrt(rNum. / rDen.)
    # phase, Kitagawa equation 21
    Y. <- (Kf. - Ku. * B. / Aw.) * TVFRG * A1
    X. <- (Ku. * B. / Aw. * TVFRG - A2) * (Kf. * TVFRG - A2) + A1 * A1
    Phs. <- atan2(-1*Y.,-1*X.)
    # params?
    # attributes?
    # message?
    toret <- cbind(omega, Amp., Phs.)
    return(toret)
  }