#' determine the number of clusters
#'
#' @param D distance matrix
#' @param epi threshold for cutting off
#' @return number of clusters
#' @examples
#' data(sim_dat)
#' dat <- extractInfo(sim_dat)
#' data_counts <- dat$count_matrix
#' bp <- BiocParallel::MulticoreParam(4)
#' D_c <- calD(data_counts,bp)
#' detK(D_c)

#' @export

detK <- function(D, epi = 1)
{
     
    tmp <- vapply(2:9,function(x) clusHelper(D,x), c(1,1))
    intra <- tmp[1,]
    inter <- tmp[2,]
    s <- intra / inter
    
    ss <- s
    
    if(min(ss) < epi){
    K <- which(ss < epi)[1] + 1
        }else{
        K <- 9
        }
    return(K)
}

#' function to get intra and inter distance for clusters
#'
#' @param D distance matrix
#' @param i number of clusters
#' @return vector of intra and inter distance
#' @export
clusHelper <- function(D,i){
    clusRes <- pam(D,i,diss = TRUE)
    intra <- as.numeric(clusRes$objective[1])
    x <- clusRes$id.med
    inter <- mean(D[x,x])
    return(c(intra,inter))
}
