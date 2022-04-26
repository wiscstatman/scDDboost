
#' calculate PDD when add random noise in distance matrix
#'
#' @param data normalized preprocessed transcripts
#' @param cd condition label
#' @param K number of subgroups
#' @param D distance matrix of cells
#' @param a shape param for weights
#' @param sz size factors 
#' @param hp hyper parameters for EBSeq
#' @param Posp parition patterns
#' @param iter max number of iterations for EM in EBSeq
#' @param REF refinement relation matrix
#' @param stp1 step size of hyperparameter alpha (shared by all units) in one step EM
#' @param stp2 step size of hyperparameter beta (unit specific) in one step EM
#' @return posterior probabilities under random distance matrix
#' @keywords internal
#' @export

pddRandom <- function(data, cd, K, D, a, sz, hp, Posp, iter, REF, stp1, stp2){

    cstar <- genRClus(D,a,K)
    gcl <- seq_len(nrow(data))
    tmp <- getZ1Z2(cstar,cd)
    z1 <- tmp[[1]]
    z2 <- tmp[[2]]
    res <- EBS(data,cstar,gcl,sz,iter,hp,Posp,stp1,stp2)
    DE <- res$DEpattern
    PDD <- pddAggregate(z1,z2,Posp,DE,K,REF)
    PDD
}
