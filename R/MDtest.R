MDtest <- 
    function(object, group, nperm=1000, pvalue.only=TRUE)
{
    if(!(is.matrix(object))) 
        stop("'object' must be a matrix where rows are features 
            and columns are samples")

    if(is.null(group)) 
        stop("'group' must be a vector indicating group association. 
            Possible values are 1 and 2")

    nv <- ncol(object)

    if(!is.logical(pvalue.only)) 
        stop("'pvalue.only' must be logical")

    if(length(group) != nv) 
        stop("length of 'group' must equal the number of columns in 'object'")

    if(sum(group %in% c(1,2)) < nv) 
        stop("all members in 'group' must have values 1 or 2")

    if((sum(group == 1) < 3) | (sum(group == 2) < 3)) 
        stop("there are less than 3 samples in at least one group")

    object <- object[,c(which(group == 1),which(group == 2))]
    nv1 <- sum(group == 1)
    objt <- aperm(object, c(2,1))
    Wmat <- as.matrix(dist(objt, method="euclidean", diag=TRUE, 
        upper=TRUE, p=2))
    gr <- graph_from_adjacency_matrix(Wmat, weighted=TRUE, mode="undirected")
    V(gr)[c(1:nv1)]$color <- "green"
    V(gr)[c((nv1+1):nv)]$color <- "red"
    MST <- mst(gr)
    KSranking <- HDP.ranking(MST)
    domain <- V(MST)$color
    D_perm <- array(0,c(1,nperm))

    for(itr in 1:nperm) 
    { 
        randperm <- sample(domain, replace=FALSE)
        di <- array(0, c(1,nv))
        for (i in 1:nv)
        {
        ri <- sum((KSranking[1:i]*(randperm[KSranking[1:i]] == "green"))^0.25) / sum((KSranking*(randperm[KSranking] == "green"))^0.25)
        si <- sum(randperm[KSranking[1:i]] == "red") / (nv-nv1)
        di[i] <- ri-si
        }
    D_perm[itr] <- sum(di)
    }

    di <- array(0, c(1,nv))

    for (i in 1:nv)
    {
        ri <- sum((KSranking[1:i]*(domain[KSranking[1:i]] == "green"))^0.25) / sum((KSranking*(domain[KSranking] == "green"))^0.25)
        si <- sum(domain[KSranking[1:i]] == "red") / (nv-nv1)
        di[i] <- ri-si
    }
    D_obs <- sum(di)
    pvalue <- (sum(abs(D_perm) >= abs(D_obs)) + 1) / (length(D_perm) + 1)
    if(pvalue.only) return(pvalue)
    if(!pvalue.only) return(list("statistic"=D_obs,"perm.stat"=D_perm,"p.value"=pvalue))
}
