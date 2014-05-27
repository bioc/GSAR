HDP.ranking <- 
    function(mst)
{
    rr <- farthest.nodes(mst, directed=FALSE, unconnected=TRUE)
    root <- floor(rr[1])
    terminal_nodes <- which(igraph::degree(mst) == 1)
    ltn <- length(terminal_nodes) - 1
    tn <- terminal_nodes
    tn <- tn[tn != root]
    sp <- get.shortest.paths(mst, root, to=tn)$vpath
    path_len <- shortest.paths(mst)
    break_ties <- path_len[root, tn] / max(path_len)
    depth <- array(0, c(1,ltn))
    KSranks <- root
    for(k in 1:ltn)    depth[k] <- length(sp[[k]])
    md <- max(depth) 
    adjusted_depth <- depth + break_ties
    col_nodes <- array(0, c(1,ltn))
    alphabets <- rep("",ltn)

    for (col in seq(1,md,by=1))
    {
        for(row in seq(1,ltn,by=1))
            col_nodes[row] <- sp[[row]][col]
    fcn <- factor(col_nodes)
    collevels <- levels(fcn)
    llev <- length(collevels)
    if (llev > 1) 
    {
        mpg <- tapply(adjusted_depth,fcn,max)
        sortmpg <- sort(mpg, decreasing=FALSE, index.return=TRUE)
        smpg <- sortmpg$ix
        sorted_levels <- collevels[smpg]
        for (lind in seq(1,length(smpg),by=1)) 
        {
            alphabets[which(col_nodes==sorted_levels[lind])] <- 
                paste(alphabets[which(col_nodes==sorted_levels[lind])], 
            letters[lind], sep="")
        }
    }
    }
    newranks <- sort(alphabets, decreasing=FALSE, index.return=TRUE)
    spm <- as.matrix(sp)
    sp_new <- spm[newranks$ix,]
    sp_new <- as.matrix(sp_new)

    for (k in 1:ltn)
    {
        len <- length(sp_new[[k]])
        for (u in 1:len)
        {
            if (sum(KSranks == sp_new[[k]][u]) == 0)
            {
                KSranks <- c(KSranks,sp_new[[k]][u])
            }
        }
    }
    KSranks
}
