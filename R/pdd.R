
#' calculate posterior probabilities of a gene to be differential distributed
#'
#' @param data normalized preprocessed transcripts
#' @param cd conditions label
#' @param ncores number of cores for parallel computing
#' @param D distance matrix of cells or cluster of cells or a given clustering
#' @param epi tol for change of validity score in determining number of clusters
#' @param random boolean indicator of whether randomzation has been been implemented on distance matrix
#' @param norm boolean indicator of whether the input expression data is normalized
#' @param Upper bound for hyper parameters optimization
#' @param nrandom number of random generated distance matrix
#' @param iter max number of iterations for EM
#' @param reltol relative tolerance for optim on weighting paramters
#' @param stp1 step size of hyperparameter alpha (shared by all units) in one step EM
#' @param stp2 step size of hyperparameter beta (unit specific) in one step EM
#' @param K number of subtypes, could be user specified or determined internally(set to 0)
#' @return posterior probabilities of a gene to be differential distributed

#' @examples
#' data(sim_dat)
#' dat = extractInfo(sim_dat)
#' data_counts = dat$count_matrix
#' cd = dat$condition
#' D_c = cal_D(data_counts,4)
#' pDD = PDD(data_counts,cd,2,D_c)
#' @export


PDD <- function(data, cd, ncores,D, random = TRUE, norm = TRUE, epi = 1, Upper = 1000, nrandom = 50, iter = 20,reltol = 1e-3, stp1 = 1e-6, stp2 = 1e-2, K = 0){
    
    G <- nrow(data)
    
    rname <- rownames(data)
    if(is.null(rname)){
        rname <- vapply(seq_len(G),function(x) paste0("gene",x),"string")
    }
    rs <- rowSums(data)
    zGene <- which(rs == 0)
    msg <- paste0(length(zGene), " genes are all zero counts, not being considered in DD analysis")
    message(msg)
    
    selected <- which(rs > 0)
    
    data <- data[selected,]
    
    gcl <- seq_len(nrow(data))
    
    
    if(norm)
    {
        if(is.matrix(D)){
            sz <- rep(1, ncol(D))
        }else if(is.numeric(D)){
            sz <- rep(1, length(D))
        }
    }else{
    sz <- tryCatch({MedianNorm(data)},error = function(e){
        message("sizeFactor calculation failed, try normalized data")
    })
    }
    alpha <- 0.4
    beta <- 2
    
    hp <- rep(beta, 1 + nrow(data))
    hp[1] <- alpha
    
    
    if(!random){
        if(is.matrix(D)){
            if(K == 0){
                K <- detK(D,epi)
            }
            msg <- paste0("estimated number of subtypes: ",K)
            message(msg)
            ccl <- pam(D, k = K, diss = TRUE)$clustering
        }
        else{
            if(is.numeric(D))
            {
                ccl <- D
                K <- max(ccl)
            }
            else
            {
               stop("input ccl must either be a dissimilarity measure or a partition for cells")
            }
            
        }
        
        
        Posp <- pat(K)[[1]]
        if(K >= 2){
            res <- EBS(data,ccl,gcl,sz,iter,hp,Posp,stp1,stp2)
            DE <- res$DEpattern
        }
        else if(K == 1){
            message("There is only one cluster, no postive")
            return(rep(0,nrow(data)))
        }
        n1 <- table(cd)[1]
        z1<-rep(0, K)
        z2<-rep(0, K)
        for(i in seq_len(K)){
            ##current index
            cur<-which(ccl==i)
            z1[i]<-length(which(cur<=n1))
            z2[i]<-length(which(cur>n1))
        }
        alpha1 <- rep(1,K)
        alpha2 <- rep(1,K)
        post <- MDD(z1, z2, Posp, alpha1, alpha2)
        np <- nrow(Posp)
        REF <- g_ref(Posp)
        modified_p <- t(REF) %*% post
        PED <- DE%*%modified_p
        
        PDD <- 1 - PED
        res <- rep(0,G)
        res[selected] <- PDD
        
        return(res)
    }
    else{
        if(K == 0){
            K <- detK(D,epi)
        }
        msg <- paste0("estimated number of subtypes: ",K)
        message(msg)
        
        if(K == 1){
            message("There is only one cluster, no postive")
            return(rep(0,nrow(data)))
        }
        
        Posp <- pat(K)[[1]]
        REF <- g_ref(Posp)
        
        # MLE for random weighting parameter
        a <- rwMLE(D,reltol)
        
        bp <- BiocParallel::MulticoreParam(ncores)
        result <- bplapply(seq_len(nrandom), function(i) {PDD_random(data, cd, K, D, a, sz, hp, Posp, iter, REF,stp1,stp2)}, BPPARAM = bp)
        
        
        boot <- matrix(0,nrow=length(result[[1]]),ncol = nrandom)
        for(i in seq_len(nrandom)){
            boot[,i] <- result[[i]]
        }
        
        PDD <- rowSums(boot) / nrandom
        res <- rep(0,G)
        res[selected] <- PDD
        
        names(res) <- rname
        
        return (res)
    }
    
}
