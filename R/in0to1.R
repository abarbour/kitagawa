#' Check if in [0,1]
#' 
# \emph{This function is not likely to be needed by the user.}
#' 
#' @name .in0to1
#' @rdname in0to1
#' @export
#' 
#' @param X   something to be checked (vector, scalar, ...)
#'
# @return NULL
#' 
#' @author Andrew Barbour <andy.barbour@@gmail.com>
#' 
#' @examples
#' \dontrun{
#' .in0to1(1:10) # error
#' .in0to1(NULL) # error
#' .in0to1(c(1:10,NULL)) # error
#' .in0to1(NA) # error
#' .in0to1(c(1:10,NA)) # error
#' }
.in0to1 <-
  function(X){stopifnot((X>=0) & (X<=1))}