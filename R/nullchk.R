nullchk <-
function(X){stopifnot(!is.null(X) & !(NA %in% X))}
