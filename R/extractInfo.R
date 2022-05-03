#' extract count matrix from SingleCellExperiment object
#'
#'@param data SingleCellExperiment object
#'@return list of count matrix and condition vector
#'
#' @examples
#' data(sim_dat) 
#' dat <- extractInfo(sim_dat)

#' @export

extractInfo <- function(data){
    if(!is(data,"SingleCellExperiment")){
        stop("input data must be SingleCellExperiment Object")
    }
    data_counts <- assays(data)$count
    cd <- colData(data)$conditions
    res <- list()
    res[['count_matrix']] <- data_counts
    res[['condition']] <- cd
    return(res)
}
