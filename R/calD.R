
#' calculate distance matrix
#'
#' @param data transcripts
#' @param bp bioc parallel parameter
#' @return distance matrix
#' @examples
#' data(sim_dat)
#' dat <- extractInfo(sim_dat)
#' data_counts <- dat$count_matrix
#' bp <- BiocParallel::MulticoreParam(4)
#' D_c <- calD(data_counts,bp)

#' @export



calD <- function(data,bp){
    nc <- ncol(data)
    nr <- nrow(data)
    
    geneMax <- apply(data,1,max)
    tmp <- which(geneMax > 0)
    cl <- gCl(data[tmp,], bp)
    
    geneK <- apply(cl,1,max)
    tmp <- which(geneK > 1)
    f_cl <- cl[tmp,]
    geneK <- geneK[tmp]
    
    
    ng <- nrow(f_cl)
    
    
    D_E <- as.matrix(dist(t(f_cl),method = "manhattan")) / nrow(f_cl)
    
    m_ <- max(D_E)
    
    D_E <- D_E / m_
    
    D_cor <- (1 - cor(data)) / 2
   
    w1 <- 1 / sd(D_E)
    
    w2 <- 1 / sd(D_cor)
    
    total <- w1 + w2
    
    w1 <- w1 / total
    
    w2 <- w2 / total
    
    D_c <- w2*D_cor + w1*D_E
    
    return(D_c)
}
