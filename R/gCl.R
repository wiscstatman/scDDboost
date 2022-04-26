#' gene_level cluster
#'
#' @param data transcripts
#' @param bp bioc parallel parameter
#' @return return a matrix whose row represent gene specific cluster
#' @export

gCl <- function(data, bp){
    nr <- nrow(data)
    Phi_mdf <- rep(1, nr)
    bt <- rep(1, nr)
    tryCatch(
    {
        sz <- rep(1, ncol(data))
        MV<-CalcMV(data,Sizes = sz, Plot = FALSE)
        Phi_mdf<-MV$Phi_mdf
        Q_mdf<-MV$Q_mdf
        bt<-1/Q_mdf-1
    },
    error = function(w) {
     message("estimation of hyper parameter failed, try naively assigned parameters")
        }
    ,finally = {
    clus <- bplapply(seq_len(nr),function(i) MCP(data[i,],1,c(Phi_mdf[i],bt[i]))+1,BPPARAM = bp)
    cl <- matrix(0, nrow = nrow(data), ncol = ncol(data))
    for(i in seq_len(length(clus))){
        cl[i,] <- clus[[i]]
    }
    return(cl)
    })
}


