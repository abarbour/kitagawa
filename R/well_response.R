#' Calculate the pressure/strain response spectrum for given formation properties
#'
#' This is the primary function which calculates the theoretical, complex
#' well response, namely \strong{Equation 17} in Kitagawa et al (2011).
#' The results, however, are expressed as amplitude and phase.
#' 
#' \strong{The response depends strongly on the physical properties
#' given.
#' Default values are assumed where resonable, mostly that the pore-fluid
#' is water, but considerable care should be invested in the choice of
#' parameters, unless the function is used in an optimization scheme.}
#' 
#' Assumed values are:
#' \describe{
#' \item{\code{Avs.}}{1}
#' \item{\code{Aw.}}{1}
#' \item{\code{rho.}}{\eqn{1000 [kg/m^3]}, the density of water}
#' \item{\code{Kf.}}{\eqn{2.2e9 [Pa]}, the bulk modulus of water}
#' \item{\code{grav.}}{\eqn{9.81 [m/s^2]}, average gravitational force on Earth}
#' }
#' 
#' Note that Skempton's coefficient, \code{B.}, is bounded inclusively
#' within \eqn{[0,1]}; an error is thrown if it's not.
#'
#' @name well_response
#' @export
#' 
#' @param omega  frequency,  (see freq.units)
#' @param T.     effective aquifer transmissivity \eqn{[m^2/s]}
#' @param S.     well storativity,  \eqn{[unitless]}
#' @param Vw.    well volume,	 \eqn{[m^3]}
#' @param Rs.    radius of screened portion,  \eqn{[m]}
#' @param Ku.    undrained bulk modulus,  \eqn{[Pa]}
#' @param B.     Skempton's coefficient,  \eqn{[unitless, bounded]}
#' @param Avs.   amplification factor for volumetric strain \eqn{E_{kk,obs}/E_{kk}},  \eqn{[]}
#' @param Aw.    amplification factor of well volume change for \eqn{E_{kk}},  \eqn{[]}
#' @param rho.   fluid density \eqn{[kg/m^3]}
#' @param Kf.    bulk modulus of fluid,  \eqn{[Pa]}
#' @param grav.  local gravitational acceleration \eqn{[m/s^2]}
#' @param freq.units  set what the units of frequency (omega) are: \code{"rad_per_sec"} (default, \code{NULL}), or \code{"Hz"}
# @param ...    additional parameters
#'
#' @return Matrix with three columns: radial frequency, amplitude, and phase 
#' [\eqn{\omega}), \eqn{A_\alpha (\omega)}, \eqn{\Phi_\alpha (\omega)}]
#' where the units of \eqn{\omega} will be radians per second,
#' \eqn{A_\alpha (\omega)} in meters per strain, 
#' and \eqn{\Phi_\alpha (\omega)} in radians.
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @references Kitagawa, Y., S. Itaba, N. Matsumoto, and N. Koisumi (2011),
#' Frequency characteristics of the response of water pressure in a closed well to volumetric strain in the high-frequency domain,
#' \emph{J. Geophys. Res.}, \strong{116}, B08301, doi:10.1029/2010JB007794.
#' 
#' @references \url{http://www.agu.org/pubs/crossref/2011/2010JB007794.shtml}
#'
#' @seealso 
#' \code{\link{sensing_volume}} to estimate the volume of water, and
#' \code{\link{kitplot}} to plot the results.
#'
#' \code{\link{open_well_response}} for modeling response of an open (exposed)
#' well.
#' 
#' @examples
#' #### dummy example
#' well_response(1:10, T.=1, S.=1, Vw.=1, Rs.=1, Ku.=1, B.=1)
#' 
#' #### a more physically realistic calculation:
#' # Physical params applicable for B084 borehole
#' # (see: http://pbo.unavco.org/station/overview/B084/ for details)
#' #
#' Rc <- 0.0508   # m, radius of water-sensing (2in)
#' Lc <- 146.9    # m, length of grouted region (482ft)
#' Rs <- 3*Rc     # m, radius of screened region (6in)
#' Ls <- 9.14     # m, length of screened region (30ft)
#' #
#' # calculate the sensing volume for the given well parameters
#' Volw <- sensing_volume(Rc, Lc, Rs, Ls) # m**3, ~= 1.8
#' #
#' Frqs <- 10**seq.int(from=-4,to=0,by=0.1) # log10-space
#' head(Rsp <- well_response(omega=Frqs, T.=1e-6, S.=1e-5, 
#' Vw.=Volw, Rs.=Rs, Ku.=40e9, B.=0.2, freq.units="Hz"))
#' #
#' kitplot(Rsp)
#'
well_response <- function(omega, T., S., Vw., Rs., Ku., B., 
           Avs.=1,
           Aw.=1,
           rho.=1000, 
           Kf.=2.2e9,
           grav.=9.81,
           freq.units=NULL) UseMethod("well_response")

#' @rdname well_response
#' @docType methods
#' @method well_response default
#' @S3method well_response default
well_response.default <- function(omega, T., S., Vw., Rs., Ku., B.,
           Avs.=1,
           Aw.=1,
           rho.=1000, 
           Kf.=2.2e9,
           grav.=9.81,
           freq.units=NULL){
    #
    # Enforce units of omega to be radians/sec
    fc <- switch(match.arg(freq.units, c("rad_per_sec","Hz")),
                 rad_per_sec=1,
                 Hz=2*pi
    )
    omega <- fc*omega
    #
    # Alpha function
    Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
    #   A1, and A2 (functions of Phi and Psi, calculated internally)
    Amat <- alpha_constants(Alpha., c.type="A")
    stopifnot(ncol(Amat)==7)
    rm(Alpha.)  # cleanup
    #  A1,2 are in Mod(A.[,6:7]) 
    A12 <- matrix(Mod(Amat[,6:7]),ncol=2)  # is complex, but imag is zero, so == abs
    stopifnot(ncol(A12)==2)
    rm(Amat)    # cleanup
    A1 <- A12[,1]
    A2 <- A12[,2]
    rm(A12)     # cleanup
    #
    rhog <- rho. * grav.
    TVFRG <- 2 * pi * T. / omega / Vw. / rhog
    #
    #check Skemptons coefficient: it should be bound to [0,1]
    .in0to1(B.)
    #
    # calculate amp and phase of response
    #
    tmpd. <- Ku. * B. / Aw. * TVFRG - A2
    rNum. <- tmpd. * tmpd. + A1 * A1
    #
    tmpd. <- Kf. * TVFRG  -  A2
    rDen. <- tmpd. * tmpd. + A1 * A1
    rm(tmpd.)
    ##
    ## complex response EQ 17
    ## cNum <- complex(real=(Ku. * B. / Aw. * TVFRG - A2), imaginary=A1)
    ## cDen <- complex(real=(Kf. * TVFRG  -  A2), imaginary=A1)
    ## cResp <- -1 * Kf. * Aw. / Avs. / rhog * cNum / cDen
    ## amplitude
    ## Amp. <- Mod(cResp)
    ## Phs. <- Arg(cResp)
    ##
    ## amplitude, Kitagawa equation 20
    ##
    Amp. <- Kf. * Aw. / Avs. / rhog * sqrt(rNum. / rDen.)
    ##
    ## phase, Kitagawa equation 21
    ##
    Y. <- (Kf. - Ku. * B. / Aw.) * TVFRG * A1
    X. <- (Ku. * B. / Aw. * TVFRG - A2) * (Kf. * TVFRG - A2) + A1 * A1
    Phs. <- atan2(-1*Y.,-1*X.)
    #
    # params?
    # attributes?
    # message?
    #
    # return results
    toret <- cbind(omega, Amp., Phs.)
    return(toret)
  }

#' Calculate the pressure/strain response spectrum for given formation properties (OPEN well)
#' 
#' Assumed values are:
#' \describe{
#' \item{\code{Avs.}}{1}
#' \item{\code{Aw.}}{1}
#' \item{\code{rho.}}{\eqn{1000 [kg/m^3]}, the density of water}
#' \item{\code{Kf.}}{\eqn{2.2e9 [Pa]}, the bulk modulus of water}
#' \item{\code{grav.}}{\eqn{9.81 [m/s^2]}, average gravitational force on Earth}
#' }
#' 
#' @name open_well_response
#' @export
#' 
#' @param omega  frequency,  (see freq.units)
#' @param T.     effective aquifer transmissivity \eqn{[m^2/s]}
#' @param S.     well storativity,  \eqn{[unitless]}
#' @param Vw.    well volume,	 \eqn{[m^3]}
#' @param Rs.    radius of screened portion,  \eqn{[m]}
#' @param Ku.    undrained bulk modulus,  \eqn{[Pa]}
#' @param B.     Skempton's coefficient,  \eqn{[unitless, bounded]}
#' @param Avs.   amplification factor for volumetric strain \eqn{E_{kk,obs}/E_{kk}},  \eqn{[]}
#' @param Aw.    amplification factor of well volume change for \eqn{E_{kk}},  \eqn{[]}
#' @param rho.   fluid density \eqn{[kg/m^3]}
#' @param Kf.    bulk modulus of fluid,  \eqn{[Pa]}
#' @param grav.  local gravitational acceleration \eqn{[m/s^2]}
#' @param freq.units  set what the units of frequency (omega) are: \code{"rad_per_sec"} (default, \code{NULL}), or \code{"Hz"}
#' @param model  character; use the model of Liu et al (1989), or Cooper et al (1965)
# @param ...    additional parameters
#'
#' @return Matrix with three columns: radial frequency, amplitude, and phase 
#' [\eqn{\omega}), \eqn{A_\alpha (\omega)}, \eqn{\Phi_\alpha (\omega)}]
#' where the units of \eqn{\omega} will be radians per second,
#' \eqn{A_\alpha (\omega)} in meters per strain, 
#' and \eqn{\Phi_\alpha (\omega)} in radians.
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @seealso \code{\link{well_response}}
#'
open_well_response <- function(omega, T., S., Vw., Rs., Ku., B., 
         Avs.=1,
         Aw.=1,
         rho.=1000, 
         Kf.=2.2e9,
         grav.=9.81,
         freq.units=NULL,
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
         freq.units=NULL,
         model=c("liu","cooper")){
  Alpha. <- omega_constants(omega, c.type="alpha", S.=S., T.=T., Rs.=Rs.)
  Kmat <- alpha_constants(Alpha., c.type="Kel")
  # alpha Kel0 Kel1, Kels are complex
  Kel0. <- Kmat[,2]
  kei. <- Im(Kel0.)
  ker. <- Re(Kel0.)
  omegsq <- omega**2
  #
  model <- match.arg(model)
  Hw. <- d. <- 1 # fix this!
  #
  if (model=="liu"){
	#
    # from Liu et al (1989), expressed in
    # Roeloffs 1998 eq 22
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
    # from Cooper et al (1965), expressed in
    # Roeloffs 1998 eq 20
    #
    #
    # the effective height of the water column in the well
    He. <- Hw. + 3*d./8 
    # d is the aquifer thickness
    # Hw is the height of the water column above the upper limit of the aquifer
    cT. <- omega * Rs.**2 / 2 / T.
    A. <- 1 - cT. * kei. + omegsq * He. / grav.
    B. <- cT. * ker.
    # the amplification of water level in the well relative to 
    # pressure head in the aquifer
	wellresp <- sqrt(A.**2 + B.**2)  # |w/h|
	# eq 10
	# h/e_kk
	# h/e_areal ?
	# wellresp <- wellresp * -2 * mu * B. * (1+nu_u) /(rho. * grav * 3 * (1-2 * nu_u)) # | w/h * h/e_kk|
	omega_peak <- sqrt(grav./He.) # eq 21
  }
  toret <- cbind(omega, wellresp)
  return(toret)
}
