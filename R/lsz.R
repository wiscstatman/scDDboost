
#' index of DD genes under FDR control
#'
#' @param pDD probability of genes being DD
#' @param FDR fdr to be controlled
#' @return index of positive genes
#' @examples
#' p_dd = c(0.01,0.99,0.7,0.5)
#' lsz(p_dd)

#' @export

lsz = function(pDD, FDR=0.01)

{
    
    ee <- 1-pDD
    
    oe <- sort(ee)
    
    or = order(ee)
    
    ff <- cumsum(oe)/seq_len(length(oe))
    
    return(or[which(ff < FDR)])
    
}
