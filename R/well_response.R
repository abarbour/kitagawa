#' Spectral response for a sealed well
#'
#' This is the primary function for a sealed well.
#' It calculates the theoretical, complex
#' well response, namely \strong{Equation 17} in Kitagawa et al (2011).
#' The results, however, are expressed as amplitude and phase.
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
#' \code{Avs}  \tab 1 \tab \tab amplification factor for volumetric strain\cr
#' \code{Aw}   \tab 1 \tab \tab amplification factor for water well\cr
#' }
#' 
#' \emph{Note that Skempton's coefficient, \code{B.}, is bounded inclusively
#' within \eqn{[0,1]}; an error is thrown if it's not.
#' }
#' 
#' Not all parameters need to given, as some properties can be assumed 
#' (say, for water).  \emph{The parameters which do not end in \code{.} do
#' not need to be specified (they may be excluded).}
#'
#' @name well_response
#' @export
#' 
#' @param omega  frequency, (see \code{freq.units})
#' @param T.     effective aquifer transmissivity \eqn{[m^2/s]}
#' @param S.     well storativity,  \eqn{[unitless]}
#' @param Vw.    well volume,	 \eqn{[m^3]}
#' @param Rs.    radius of screened portion,  \eqn{[m]}
#' @param Ku.    undrained bulk modulus,  \eqn{[Pa]}
#' @param B.     Skempton's coefficient,  \eqn{[unitless, bounded]}
#' @param Avs   amplification factor for volumetric strain \eqn{E_{kk,obs}/E_{kk}},  \eqn{[]}
#' @param Aw    amplification factor of well volume change for \eqn{E_{kk}},  \eqn{[]}
#' @param rho   fluid density \eqn{[kg/m^3]}
#' @param Kf    bulk modulus of fluid,  \eqn{[Pa]}
#' @param grav  local gravitational acceleration \eqn{[m/s^2]}
#' @param freq.units  set the units of \code{omega}
#'
#' @return Matrix with three columns: radial frequency, amplitude, and phase 
#' [\eqn{\omega}), \eqn{A_\alpha (\omega)}, \eqn{\Phi_\alpha (\omega)}]
#' where the units of \eqn{\omega} will be as they were input,
#' \eqn{A_\alpha (\omega)} in meters per strain, 
#' and \eqn{\Phi_\alpha (\omega)} in radians.
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#'
#' @seealso 
#' \code{\link{sensing_volume}} to estimate the volume \code{Vw.}
#' 
#' \code{\link{kitplot}} to plot the results.
#'
#' \code{\link{open_well_response}} for modeling an open well.
#' 
#' \code{\link{constants}} for assumed constants.
#' 
#' @family WellResponseFunctions
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
           Avs,
           Aw,
           rho, 
           Kf,
           grav,
           freq.units=c("rad_per_sec","Hz")) UseMethod("well_response")

#' @rdname well_response
#' @method well_response default
#' @S3method well_response default
well_response.default <- function(omega, T., S., Vw., Rs., Ku., B.,
           Avs,
           Aw,
           rho, 
           Kf,
           grav,
           freq.units=c("rad_per_sec","Hz")){
    # Enforce units of omega to be radians/sec
    freq.units <- match.arg(freq.units)
    fc <- switch(freq.units, rad_per_sec=1, Hz=2*pi)
    omega <- fc*omega
    #
    # Setup constants
    defA <- 1
    if (missing(Avs)) Avs <- defA
    if (missing(Aw)) Aw <- defA
    const <- kitagawa::constants(FALSE)
    if (missing(rho)) rho <- const$water$density
    if (missing(Kf)) Kf <- const$water$bulkmod
    if (missing(grav)) grav <- const$gravity
    rhog <- rho*grav
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
    TVFRG <- 2 * pi * T. / omega / Vw. / rhog
    #
    #check Skemptons coefficient: it should be bound to [0,1]
    .in0to1(B.)
    #
    # calculate amp and phase of response
    #
    tmpd. <- Ku. * B. / Aw * TVFRG - A2
    rNum. <- tmpd. * tmpd. + A1 * A1
    #
    tmpd. <- Kf * TVFRG  -  A2
    rDen. <- tmpd. * tmpd. + A1 * A1
    rm(tmpd.)
    ##
    ## complex response EQ 17
    ## cNum <- complex(real=(Ku. * B. / Aw * TVFRG - A2), imaginary=A1)
    ## cDen <- complex(real=(Kf * TVFRG  -  A2), imaginary=A1)
    ## cResp <- -1 * Kf * Aw / Avs / rhog * cNum / cDen
    ## amplitude
    ## Amp. <- Mod(cResp)
    ## Phs. <- Arg(cResp)
    ##
    ## amplitude, Kitagawa equation 20
    ##
    Amp. <- Kf * Aw / Avs / rhog * sqrt(rNum. / rDen.)
    ##
    ## phase, Kitagawa equation 21
    ##
    Y. <- (Kf - Ku. * B. / Aw) * TVFRG * A1
    X. <- (Ku. * B. / Aw * TVFRG - A2) * (Kf * TVFRG - A2) + A1 * A1
    Phs. <- atan2(-1*Y.,-1*X.)
    #
    # TODO: Return complex only
    #
    # return results
    omega <- omega/fc
    toret <- list(Aquifer=list(Transmiss=T., Storativ=S., Diffusiv=T./S.),
                  Well=list(Volume=Vw., ScreenRad=Rs.),
                  Solid=list(BulkModU=Ku.,Skemp=B.),
                  Fluid=list(BulkMod=Kf, Density=rho),
                  Amplification=list(Ekk=Avs, Well=Aw),
                  Omega=list(Units=freq.units),
                  Gravity=grav,
                  Response=cbind(omega=omega, Amp.=Amp., Phs.=Phs.))
    class(toret) <- "wrsp"
    return(toret)
}

#' @rdname well_response
#' @aliases print.wrsp
#' @method print wrsp
#' @S3method print wrsp
print.wrsp <- function(X, ...){
  stopifnot(is.wrsp(X))
  message("Sealed well response:")
  Xm <- as.data.frame(X$Response)
  print(head(Xm,3))
  message("\t...")
  print(tail(Xm,3))
}
# 
# #' @rdname well_response
# #' @aliases as.data.frame.wrsp
# #' @method as.data.frame wrsp
# #' @S3method as.data.frame wrsp
# as.data.frame.wrsp <- function(X, ...){
#   df <- as.data.frame.numeric(X)
#   dfn <- as.vector(attributes(X)$dimnames[[2]])
#   print(colnames(df)) # <- dfn
#   return(df)
# }
# #' @rdname well_response
# #' @aliases data.frame.wrsp
# #' @method data.frame wrsp
# #' @S3method data.frame wrsp
# data.frame.wrsp <- as.data.frame.wrsp
