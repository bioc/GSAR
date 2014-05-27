RKStest <- 
    function(object, group, nperm=1000)
{
    if(!(is.matrix(object))) 
        stop("'object' must be a matrix where rows are features and 
            columns are samples")

    if(is.null(group)) 
        stop("'group' must be a vector indicating group association. 
            Possible values are 1 and 2")

    nv <- ncol(object)

    if(length(group) != nv) 
        stop("length of 'group' must equal the number of columns in 'object'")

    if(sum(group %in% c(1,2)) < nv) 
        stop("all members in 'group' must have values 1 or 2")

    if((sum(group == 1) < 3) || (sum(group == 2) < 3)) 
        stop("there are less than 3 samples in at least one group")

    object <- object[,c(which(group == 1),which(group == 2))]
    nv1 <- sum(group == 1)
    objt <- aperm(object, c(2,1))
    Wmat <- as.matrix(dist(objt, method="euclidean", diag=TRUE, 
        upper=TRUE, p=2))
    gr <- graph.adjacency(Wmat, weighted=TRUE, mode="undirected")
    V(gr)[c(1:nv1)]$color <- "green"
    V(gr)[c((nv1+1):nv)]$color <- "red"
    mst <- minimum.spanning.tree(gr)
    sp <- apply(shortest.paths(mst), 1, max)
    radius <- min(sp)
    center <- which(sp == radius)
    if (length(center) > 1) center <- center[1]
    ranktree <- sort(shortest.paths(mst)[, center], decreasing=FALSE, 
        index.return=TRUE)
    radial_ranking <- ranktree$ix
    domain <- V(mst)$color
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
    pvalue <- (sum(D_perm > D_obs) + 1) / (length(D_perm) + 1)
    list("statistic"=D_obs,"perm.stat"=D_perm,"p.value"=pvalue)
}
