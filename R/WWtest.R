WWtest <- 
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

    object <- object[,c(which(group == 1), which(group == 2))]
    nv1 <- sum(group == 1)
    objt <- aperm(object, c(2,1))
    Wmat <- as.matrix(dist(objt, method="euclidean", diag=TRUE, 
        upper=TRUE, p=2))
    gr <- graph.adjacency(Wmat, weighted=TRUE, mode="undirected")
    V(gr)[c(1:nv1)]$color <- "green"
    V(gr)[c((nv1+1):nv)]$color <- "red"
    mst <- minimum.spanning.tree(gr)
    domain <- V(mst)$color
    runs <- array(0,c(1,nperm))

    for(itr in 1:nperm) 
    { 
        randperm <- sample(domain, replace = FALSE)
        mst2 <- mst
        V(mst2)$color <- randperm
        mstWM <- get.adjacency(mst, type="lower", attr="weight", sparse=FALSE)
        edgeind <- which(mstWM != 0, arr.ind = TRUE, useNames = FALSE) 
        runs[itr] <- 1 + 
            sum(V(mst2)[edgeind[,1]]$color != V(mst2)[edgeind[,2]]$color)
    }

    sd_runs <- apply(runs, 1, sd)
    W_perm <- (runs - mean(runs)) / sd_runs
    mstWM <- get.adjacency(mst, type = "lower", attr="weight", sparse=FALSE)
    edgeind <- which(mstWM != 0, arr.ind=TRUE, useNames=FALSE) 
    runs_obs <- 1 + sum(V(mst)[edgeind[,1]]$color != V(mst)[edgeind[,2]]$color)
    W_obs <- (runs_obs - mean(runs)) / sd_runs 
    pvalue <- (sum(W_perm < W_obs) + 1) / (length(W_perm) + 1)
    list("statistic"=W_obs,"perm.stat"=W_perm,"p.value"=pvalue)
}
