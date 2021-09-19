#' extract count matrix from SingleCellExperiment object
#'
#'@param data SingleCellExperiment object
#'@return count matrix
#'


scToMatrix = function(data){
    if(class(data)[1] != 'SingleCellExperiment'){
        stop("input data must be SingleCellExperiment Object")
    }
    data_counts = assays(data)$count
    return(data_counts)
}
