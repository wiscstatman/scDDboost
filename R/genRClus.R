
#' generate random clusterings
#'
#' @param D distance matrix of cells
#' @param a paramter for weights
#' @param K number of subtypes
#' @return random generated clustering of cells

genRClus = function(D,a,K){
    n = ncol(D)
    e <- rgamma(n,shape= a / 2, rate= a )
    bar = D/outer(e,e,"+")
    #dst.star <- as.dist(bar)
    cstar = pam(bar, k = K, diss = TRUE)$clustering
    return(cstar)
}


