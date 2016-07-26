RKStest <- 
    function(object, group, mst.order=1, nperm=1000, pvalue.only=TRUE)
{
    if(!(is.matrix(object))) 
        stop("'object' must be a matrix where rows are features and 
            columns are samples")

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

    if((sum(group == 1) < 3) || (sum(group == 2) < 3)) 
        stop("there are less than 3 samples in at least one group")

    if(mst.order > 5)
    {
        warning("'mst.order' cannot be greater than 5. Larger values are 
            reduced to 5")
        mst.order <- 5
    }

    object <- object[,c(which(group == 1),which(group == 2))]
    nv1 <- sum(group == 1)
    objt <- aperm(object, c(2,1))
    Wmat <- as.matrix(dist(objt, method="euclidean", diag=TRUE, 
        upper=TRUE, p=2))
    gr <- graph_from_adjacency_matrix(Wmat, weighted=TRUE, mode="undirected")
    MST <- mst(gr)
    if(mst.order>1) 
    {
       for(k in 1:(mst.order-1)) 
       {
           mst.matrix <- as_adjacency_matrix(MST, attr="weight", sparse=FALSE)
           distmat2 <- Wmat - mst.matrix
           gr2 <- graph_from_adjacency_matrix(distmat2, weighted=TRUE, mode="undirected")
           second.mst <- mst(gr2)
           mst2.matrix <- as_adjacency_matrix(second.mst, attr="weight", sparse=FALSE)
           mst.matrix <- mst.matrix + mst2.matrix
           MST <- graph_from_adjacency_matrix(mst.matrix, weighted=TRUE, mode="undirected")
       }
    }
    mst.matrix.unweighted <- as_adjacency_matrix(MST, attr=NULL, sparse=FALSE)
    mst.unweighted <- graph_from_adjacency_matrix(mst.matrix.unweighted, weighted=NULL,
        mode="undirected")
    V(MST)[c(1:nv1)]$color <- "green"
    V(MST)[c((nv1+1):nv)]$color <- "red"
    sp <- apply(shortest.paths(mst.unweighted), 1, max)
    spw <- apply(shortest.paths(MST), 1, max)
    new.sp <- sp + spw/max(spw)
    radius <- min(new.sp)
    center <- which(new.sp == radius)
    if (length(center) > 1) center <- center[1]
    ranktree <- sort(shortest.paths(mst.unweighted)[, center] + 
        (shortest.paths(MST)[, center]/max(shortest.paths(MST)[, center])),
        decreasing=FALSE, index.return=TRUE)
    radial_ranking <- ranktree$ix
    domain <- V(MST)$color
    D_perm <- array(0,c(1,nperm))

    for(itr in 1:nperm) 
    { 
        randperm <- sample(domain, replace=FALSE)
        ri <- 0
        si <- 0
        di <- array(0, c(1,nv))
        for (i in 1:nv)
        {
            ri <- sum(randperm[radial_ranking[1:i]] == "green")
            si <- sum(randperm[radial_ranking[1:i]] == "red")
            di[i] <- (ri/nv1) - (si/(nv-nv1))
        }
    D_perm[itr] <- sqrt((nv1 * (nv-nv1)) / (nv1 + (nv-nv1))) * max(abs(di))
    }

    ri <- 0
    si <- 0
    di <- array(0, c(1,nv))

    for (i in 1:nv)
    {
        ri <- sum(domain[radial_ranking[1:i]] == "green")
        si <- sum(domain[radial_ranking[1:i]] == "red")
        di[i] <- (ri/nv1) - (si/(nv-nv1))
    }
    D_obs <- sqrt((nv1 * (nv-nv1)) / (nv1 + (nv-nv1))) * max(abs(di))
    pvalue <- (sum(D_perm >= D_obs) + 1) / (length(D_perm) + 1)    
    if(pvalue.only) return(pvalue)
    if(!pvalue.only) return(list("statistic"=D_obs,"perm.stat"=D_perm,"p.value"=pvalue))
}
