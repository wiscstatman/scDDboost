#' generate reference matrix
#'
#' @param Posp possible partition of data
#' @return return a matrix indicate the refinement relation between different partitions.
#' @export


gRef <- function(Posp){
    .Call("g_ref",Posp = Posp)
}


