##
##
##
#
## check if something is NULL
nullchk <- function(X){stopifnot(!is.null(X) & !(NA %in% X))}
#
