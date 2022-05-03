
#' number of DD genes under FDR control
#'
#' @param pDD estimated probability of being DD
#' @param FDR fdr to be controlled
#' @return number of positive genes
#' @examples
#' p_dd <- c(0.1,0.99,1,0.05,0.05)
#' getSizeofDD(p_dd)
#' @export

getSizeofDD <- function(pDD, FDR=0.01)
{
    
    ee <- 1-pDD
    
    oe <- sort(ee)
    
    ff <- cumsum(oe)/seq_len(length(oe))
    
    return( sum( ff <= FDR ) )
    
}
