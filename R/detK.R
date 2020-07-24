#' determine the number of clusters
#'
#' @param D distance matrix
#' @param epi threshold for cutting off
#' @return number of clusters
#' @examples
#' data(sim_dat)
#' data_counts = assays(sim_dat)$count
#' D_c = cal_D(data_counts,4)
#' detK(D_c)

#' @export

detK = function(D, epi = 1)
{
     intra = rep(0,8)
     inter = rep(0,8)
    
     for(i in 2:9){
         clusRes = pam(D,i,diss = T)
         intra[i - 1] = as.numeric(clusRes$objective[1])
         x = clusRes$id.med
         #inter[i - 1] = sum(D[x,x]) / (i * (i - 1))
         inter[i - 1] = mean(D[x,x])
     }
    
    s = intra / inter
    
    #mins = min(s)
    #ss = s/mins - 1
    
    ss = s
    
    if(min(ss) < epi){
    K = which(ss < epi)[1] + 1
        }else{
        K = 9
        }
    return(K)
}
