##
##
##
# well volume calculation
sensing_volume <- function(rad_case, len_case, rad_screen, len_screen){
  Rc. <- rad_case
  Lc. <- len_case
  Rs. <- rad_screen
  Ls. <- len_screen
  Vw. <- Rc. * Rc. * (Lc. - Ls.) + Rs. * Rs. * Ls.
  return(pi*Vw.)
}
#
## check if something is NULL
nullchk <- function(X){stopifnot(!is.null(X) & !(NA %in% X))}
#
