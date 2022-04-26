#' determine the number of clusters
#'
#' @param D distance matrix
#' @param epi threshold for cutting off
#' @return number of clusters
#' @examples
#' data(sim_dat)
#' dat = extractInfo(sim_dat)
#' data_counts = dat$count_matrix
#' D_c = cal_D(data_counts,4)
#' detK(D_c)

#' @export

detK <- function(D, epi = 1)
{
     intra <- rep(0,8)
     inter <- rep(0,8)
    
     for(i in 2:9){
         clusRes <- pam(D,i,diss = TRUE)
         intra[i - 1] <- as.numeric(clusRes$objective[1])
         x <- clusRes$id.med
         inter[i - 1] <- mean(D[x,x])
     }
    
    s <- intra / inter
    
    ss <- s
    
    if(min(ss) < epi){
    K <- which(ss < epi)[1] + 1
        }else{
        K <- 9
        }
    return(K)
}
