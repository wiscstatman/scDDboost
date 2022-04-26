
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
    n1 <- table(cd)[1]
    n <- length(cstar)
    tmp_z1<-table(cstar[1:n1])
    tmp_z2<-table(cstar[(n1+1):n])
    z1 <- rep(0,K)
    z2 <- rep(0,K)
    z1[as.numeric(names(tmp_z1))] = as.numeric(tmp_z1)
    z2[as.numeric(names(tmp_z2))] = as.numeric(tmp_z2)
    
    alpha1 <- rep(1,K)
    alpha2 <- rep(1,K)
    post <- mdd(z1, z2, Posp, alpha1, alpha2)
    np <- nrow(Posp)
    modified_p <- t(REF) %*% post
    
    if(K >= 2){
        res <- EBS(data,cstar,gcl,sz,iter,hp,Posp,stp1,stp2)
        DE <- res$DEpattern
    }
    PED <- DE%*%modified_p
    
    
    PDD <- 1 - PED
    return(PDD)
}
