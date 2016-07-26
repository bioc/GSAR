HDP.ranking <- 
    function(object)
{
    if(!is_igraph(object)) 
        stop("'object' must be of class igraph. See package igraph for details")

    rr <- farthest.nodes(object, directed=FALSE, unconnected=TRUE)
    if(is.list(rr)) root <- floor(as.numeric(rr$vertices[1]))
    if(is.numeric(rr)) root <- floor(rr[1])
    terminal_nodes <- which(igraph::degree(object) == 1)
    ltn <- length(terminal_nodes) - 1
    tn <- terminal_nodes
    tn <- tn[tn != root]
    if(is.list(rr)) 
    sp <- all_shortest_paths(object, from=root, to=tn, mode="all")$res
    if(is.numeric(rr))
    sp <- get.shortest.paths(object, from=root, to=tn, mode="all")$vpath
    path_len <- shortest.paths(object)
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
            if(length(sp[[row]])>=col) 
			col_nodes[row] <- sp[[row]][col] 
				else col_nodes[row] <- 0
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
    as.numeric(KSranks)
}
